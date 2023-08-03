#! /usr/bin/python3


import re
import sys
import subprocess


class augustus_res():
    def __init__(self):
        self.gene = ""
        self.gff = []
        self.codingseq = ""
        self.aaseq = ""
        self.percent_support = 0.0
        self.support_exons = 0
        self.support_introns = 0
        self.five_utr = 0
        self.three_utr = 0
        self.fully = 0
        self.incompatible = 0


def run_augustus(dbf: str, hintsf: str, species="arabidopsis",
                 coding=True) -> str:
    cmd = ["augustus", "--species=" + species,
           "--hintsfile=" + hintsf, dbf]
    if coding:
        cmd += ["--codingseq=on"]
    sys.stderr.write(subprocess.list2cmdline(cmd)+"\n\n")
    p = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    out = dbf + ".hints.augustus.out"
    with open(out, "w") as outf:
        for line in p.stdout.splitlines():
            outf.write(line+"\n")
    return out


def parse_augustus(outf: str) -> list[augustus_res]:
    out = []
    with open(outf, "r") as res:
        gene = False
        for line in res:
            if gene:
                if not line.startswith("#"):  # gff
                    newres.gff.append(line.strip())
                if cds:
                    
            if line.startswith("# start gene"):
                gene = True
                newres = augustus_res()
                line = re.sub("# start gene ", "", line).strip()
                newres.gene = line
            

            
    return out
