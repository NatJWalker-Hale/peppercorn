#! /usr/bin/python3

import os
import sys
import argparse
import subprocess


class tblastnhit():
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
            if h.evalue < e and h.bitscore > bitscore:
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
                j.sstartorder <= hit.sstartorder and
                j.sendorder >= hit.sendorder  # fully subsumed
            ):
                if j.evalue < hit.evalue and j.bitscore > hit.bitscore:
                    out.remove(hit)
                else:
                    out.remove(j)
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
