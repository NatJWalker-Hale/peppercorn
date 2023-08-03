#! /usr/bin/python3


import sys
import argparse
import subprocess


def run_exonerate(dbf, queryf, model="protein2genome", trim=15):
    sys.stderr.write("Aligning proteins to genome with Exonerate\n\n")
    cmd = ["exonerate", "-m", model, "--showtargetgff", "TRUE", queryf,
           dbf]
    sys.stderr.write(subprocess.list2cmdline(cmd)+"\n\n")
    p = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    out = ".".join([dbf, queryf, "exonerate.out"])
    with open(out, "w") as outf:
        for line in p.stdout.splitlines():
            outf.write(line.strip()+"\n")
    outgff = ".".join([dbf, queryf, "exonerate.hints"])
    with open(outgff, "w") as outf:
        going = False
        result = False
        for line in p.stdout.splitlines():
            if line == "# --- START OF GFF DUMP ---":
                going = True
                result = True
            elif line == "# --- END OF GFF DUMP ---":
                going = False
            if going:
                if not line.startswith("#"):
                    rec = line.split("\t")
                    if rec[2] == "intron":
                        outf.write("\t".join([rec[0], "xnt2h", "intron",
                                              rec[3], rec[4], rec[5], rec[6],
                                              rec[7],
                                              "src=M;grep=" + queryf +
                                              ";pri=4\n"
                                              ]))
                    elif rec[2] == "cds":
                        outf.write("\t".join([rec[0], "xnt2h", "CDSpart",
                                              str(int(rec[3]) + trim),
                                              str(int(rec[4]) - trim),
                                              rec[5], rec[6],
                                              rec[7],
                                              "src=M;grep=" + queryf +
                                              ";pri=4\n"
                                              ]))
        return out, outgff, result


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("database", help="FASTA-formatted database to align \
                        queries with exonerate")
    parser.add_argument("query", help="FASTA-formatted protein queries to \
                        align")
    args = parser.parse_args()

    xnt_out, xnt_hint, _ = run_exonerate(args.database, args.query)
    print(xnt_out, xnt_hint)
