#! /usr/bin/python3


import re
import sys
import subprocess


class AugustusRes():
    """
    simple class to hold Augustus results
    """
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
        """
        write FASTA of self.codingseq (default) or self.aaseq
        """
        if coding:
            out = f">{self.gene}\n{self.codingseq}\n"
        else:
            out = f">{self.gene}\n{self.aaseq}\n"
        return out


def run_augustus(dbf: str, hintsf: str, species="arabidopsis",
                 coding=True) -> str:
    """
    run augustus using arabidopsis model and returning coding sequences by default
    """
    cmd = ["augustus", "--species=" + species,
           "--hintsfile=" + hintsf, dbf]
    if coding:
        cmd += ["--codingseq=on"]
    sys.stderr.write(subprocess.list2cmdline(cmd)+"\n\n")
    p = subprocess.run(cmd, shell=False, capture_output=True, text=True, check=True)
    out = dbf + ".hints.augustus.out"
    with open(out, "w", encoding="utf-8") as outf:
        for line in p.stdout.splitlines():
            outf.write(line+"\n")
    return out


def parse_augustus(outf: str) -> list[AugustusRes]:
    """
    parse augustus output gff
    """
    out = []
    with open(outf, "r", encoding="utf-8") as res:
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
                newres = AugustusRes()
                newres.gene = line.split(" ")[2]
    return out


def rename(results: list[AugustusRes]):
    """
    rename a list of Augustus annotations to consecutive integer format
    """
    i = 1
    for r in results:
        r.gene = f"g{i}"
        i += 1
