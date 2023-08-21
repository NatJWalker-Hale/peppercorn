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
        self.unsupport_exons = 0
        self.support_introns = 0
        self.unsupport_introns = 0
        self.fully = 0
        self.incompatible = 0

    def write_fasta(self, coding=True):
        if coding:
            out = ">%s\n%s\n" % (self.gene, self.codingseq)
        else:
            out = ">%s\n%s\n" % (self.gene, self.aaseq)
        return out


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
        cds = False
        aa = False
        for line in res:
            if line.startswith("#"):
                line = line.lstrip("# ").strip()
            else:  # gff
                newres.gff.append(line.strip().split())  # order is meaningful
            if gene:
                if line.startswith("end gene"):
                    gene = False
                    out.append(newres)
                if cds:
                    if line.endswith("]"):
                        cds = False
                        line = re.sub("\\]",
                                      "",
                                      line)
                        newres.codingseq += line
                        newres.codingseq = newres.codingseq.upper()
                    else:
                        newres.codingseq += line.strip()
                if line.startswith("coding sequence"):
                    cds = True
                    line = re.sub("coding sequence = \\[",
                                  "",
                                  line)
                    newres.codingseq += line
                if aa:
                    if line.endswith("]"):
                        aa = False
                        line = re.sub("\\]",
                                      "",
                                      line)
                        newres.aaseq += line
                    else:
                        newres.aaseq += line
                if line.startswith("protein sequence"):
                    aa = True
                    line = re.sub("protein sequence = \\[",
                                  "",
                                  line)
                    newres.aaseq += line
                if line.startswith("%"):
                    newres.percent_support = float(line.split(":")[1].strip())
                if line.startswith("CDS exons"):
                    num = int(line.split(":")[1].strip().split("/")[0])
                    denom = int(line.split(":")[1].strip().split("/")[1])
                    newres.support_exons = num
                    newres.unsupport_exons = denom
                if line.startswith("CDS introns"):
                    num = int(line.split(":")[1].strip().split("/")[0])
                    denom = int(line.split(":")[1].strip().split("/")[1])
                    newres.support_introns = num
                    newres.unsupport_introns = denom
                if line.startswith("hint groups"):
                    num = int(line.split(":")[1].strip())
                    newres.fully = num
                if line.startswith("incompatible"):
                    num = int(line.split(":")[1].strip())
                    newres.incompatible = num
            if line.startswith("start gene"):
                gene = True
                newres = augustus_res()
                newres.gene = line.split(" ")[2]
    return out


def rename(results: list[augustus_res]):
    i = 1
    for r in results:
        r.gene = "g%s" % i
        i += 1
