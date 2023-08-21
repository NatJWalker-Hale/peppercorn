#! /usr/bin/python3

import os
import sys
import argparse
import subprocess
from utils import tblastnhit, genome_range


def run_makeblastdb(dbf: str, dbtype="nucl"):
    sys.stderr.write("Making blast db\n\n")
    cmd = ["makeblastdb", "-in", dbf, "-dbtype", dbtype, "-parse_seqids"]
    sys.stderr.write(subprocess.list2cmdline(cmd) + "\n\n")
    subprocess.run(cmd, shell=False)


def run_tblastn(queryf: str, dbf: str, threads=1) -> str:
    sys.stderr.write("Searching for matches to quer(y/ies)\n\n")
    outf = queryf + ".tblastn.outfmt7"
    cmd = ["tblastn", "-query", queryf, "-db", dbf, "-outfmt", "7", "-out",
           outf, "-num_threads", str(threads)]
    sys.stderr.write(subprocess.list2cmdline(cmd) + "\n\n")
    subprocess.run(cmd, shell=False)
    return outf


def parse_tblastn(outf: str, e=1e-5, bitscore=30) -> list[tblastnhit]:
    sys.stderr.write("Parsing tblastn hits\n\n")
    out = []
    with open(outf, "r") as raw:
        for line in raw.readlines():
            if line.startswith("#"):
                continue
            line = line.strip().split()
            h = tblastnhit(input=line)
            if h.evalue < e and h.bitscore > bitscore:  # filter
                out.append(h)
    return out  # list of lists


def unique_hits(tblastn_out: list[tblastnhit]) -> list[tblastnhit]:
    out = []
    for h in tblastn_out:
        out.append(h)
    for hit in out:
        for j in out:
            if hit == j:
                continue  # exclude self
            if (
                j.sstartorder == hit.sstartorder and
                j.sendorder == hit.sendorder  # same subject span
            ):
                if j.evalue < hit.evalue and j.bitscore > hit.bitscore:
                    out.remove(hit)  # other hit is better
                else:
                    out.remove(j)
            elif (
                j.sstartorder < hit.sstartorder and
                j.sendorder > hit.sendorder  # fully subsumed
            ):
                out.remove(hit)
    return out


def join_hits(tblastn_out: list[tblastnhit],
              length=3000, over=500) -> list[tblastnhit]:
    """
    connects hits (assumed separate exons) into ranges if distance between
    hits is less than length. Therefore, length should be longer than expected
    max intron length (and if unsure, longer is better).

    This has the pleasing benefit of merging overlaps, because diff is -ve
    which is less than length.
    """
    out = []
    subjects = set([x.subject for x in tblastn_out])
    for s in subjects:
        ordered = sorted([x for x in tblastn_out if x.subject == s],
                         key=lambda x: x.sstartorder)
        or_iter = iter(ordered)
        last_hit = next(or_iter)
        for next_hit in or_iter:
            gap = next_hit.sstartorder - last_hit.sendorder
            if gap < length:
                last_hit.sendorder = next_hit.sendorder
                if gap > 0:  # not overlapping
                    if last_hit.query != next_hit.query:
                        if next_hit.query not in last_hit.queries:
                            last_hit.queries.append(next_hit.query)
                elif -gap > over:   # overlapping, want to merge
                    if last_hit.query != next_hit.query:
                        if (
                            last_hit.evalue < next_hit.evalue or
                            last_hit.bitscore > next_hit.bitscore
                        ):
                            sys.stderr.write("Replacing %s with %s\n\n" %
                                             (next_hit.query, last_hit.query))
                            next_hit.query = last_hit.query  # replace with
            else:
                out.append(last_hit)
                last_hit = next_hit
        out.append(last_hit)
    return out


def hits_to_range(tblastn_out: list[tblastnhit]) -> list[genome_range]:
    out = []
    for h in tblastn_out:
        newrange = genome_range([h.query, h.subject,
                                 h.sstartorder, h.sendorder])
        if len(h.queries) != 0:
            newrange.queries = h.queries
        out.append(newrange)
    return out


def pad_range(ranges: list[genome_range], length=3000):
    """
    pad a length on either end of a range
    """
    for r in ranges:
        r.start -= length
        r.end += length


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("database", help="FASTA-formatted nucleotide \
                        database to search")
    parser.add_argument("query", help="FASTA-formatted protein query")
    args = parser.parse_args()

    blastdbsuf = [".ndb", ".nhr", ".nin", ".nog", ".nos", ".not",
                  ".nsq", ".ntf", ".nto"]
    for s in blastdbsuf:
        if not os.path.isfile(args.database + s):
            run_makeblastdb(args.database)
    tblastn_out = run_tblastn(args.query, args.database)
    out = parse_tblastn(tblastn_out)
    print([(x.query, x.subject,
            x.sstartorder, x.sendorder) for x in unique_hits(out)])
    print([(x.query, x.subject,
            x.sstartorder, x.sendorder,
            x.queries) for x in join_hits(unique_hits(out))])
    print([(x.query, x.subject,
           x.start, x.end,
           x.queries) for x in hits_to_range(join_hits(unique_hits(out)))])
