#! /usr/bin/python3

import os
import sys
import argparse
import subprocess
from itertools import pairwise


class tblastnhit():
    """
    simple referencable data class for tblastn outfmt 7
    """
    def __init__(self, input=[]):
        self.query = input[0]
        self.subject = input[1]
        self.id = float(input[2])
        self.len = int(input[3])
        self.mm = int(input[4])
        self.gapopens = int(input[5])
        self.qstart = int(input[6])
        self.qend = int(input[7])
        self.sstart = int(input[8])
        self.send = int(input[9])
        self.evalue = float(input[10])
        self.bitscore = float(input[11])
        self.sstartorder = 0
        self.sendorder = 0
        if self.sstart > self.send:  # reverse
            self.sstartorder = self.send
            self.sendorder = self.sstart
        else:
            self.sstartorder = self.sstart
            self.sendorder = self.send


def run_makeblastdb(dbf: str, dbtype="nucl"):
    sys.stderr.write("Making blast db\n\n")
    cmd = ["makeblastdb", "-in", dbf, "-dbtype", dbtype, "-parse_seqids"]
    sys.stderr.write(subprocess.list2cmdline(cmd) + "\n\n")
    subprocess.run(cmd, shell=False)


def run_tblastn(queryf: str, dbf: str) -> str:
    sys.stderr.write("Searching for matches to quer(y/ies)\n\n")
    outf = queryf + ".tblastn.outfmt7"
    cmd = ["tblastn", "-query", queryf, "-db", dbf, "-outfmt", "7", "-out",
           outf]
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
                j.sstartorder <= hit.sstartorder and
                j.sendorder >= hit.sendorder  # fully subsumed
            ):
                out.remove(hit)
    return out


def join_hits(tblastn_out: list[tblastnhit], length=2000) -> list[tuple]:
    """
    connects hits (assumed separate exons) into ranges if distance between
    hits is less than length. Therefore, length should be longer than expected
    max intron length
    """
    out = []
    subjects = set([x.subject for x in tblastn_out])
    for s in subjects:
        ordered = sorted([x for x in tblastn_out if x.subject == s],
                         key=lambda x: x.sstartorder)
        # print([(x.sstartorder, x.sendorder) for x in ordered])
        for x in pairwise(ordered):
            if x[1].sstartorder - x[0].sendorder < length:
                if x[0].query != x[1].query:
                    if (
                        x[0].evalue < x[1].evalue or
                        x[0].bitscore > x[1].bitscore
                    ):
                        betterhit = x[0]
                    else:
                        betterhit = x[1]
                    out.append((betterhit.query, x[0].subject,
                                x[0].sstartorder, x[1].sendorder))
                else:
                    out.append((x[0].query, x[0].subject,
                                x[0].sstartorder, x[1].sendorder))
    return out


def pad_range(tblastn_out: list[tblastnhit], factor=2) -> list[tuple]:
    out = []
    subjects = set([x.subject for x in tblastn_out])
    for s in subjects:
        lens = [x.sendorder - x.sstartorder for x in tblastn_out if
                x.subject == s]
        mean_len = sum(lens) / len(lens)
        for h in tblastn_out:
            if h.subject == s:
                out.append((h.query, h.subject,
                            round(h.sstartorder - factor * mean_len),
                            round(h.sendorder + factor * mean_len)))
    return out


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
    print(join_hits(unique_hits(out)))
