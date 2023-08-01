#! /usr/bin/python3


import os
import sys
import utils
import argparse
import tblastn
import subseq
import exonerate
import augustus


def split_query(queryf):
    queryseqs = dict([x for x in utils.parse_fasta(queryf)])
    for k, v in queryseqs.items():
        with open(k + ".pep.fa", "w") as outf:
            outf.write(">%s\n" % k)
            outf.write("%s\n" % v)


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("genome", help="FASTA-formatted genomic contigs to \
                        annotate homologues")
    parser.add_argument("queries", help="FASTA-formatted protein queries for \
                        homology-based annotation")
    args = parser.parse_args()

    # run tblastn
    blastdbsuf = [".ndb", ".nhr", ".nin", ".nog", ".nos", ".not",
                  ".nsq", ".ntf", ".nto"]
    for s in blastdbsuf:
        if not os.path.isfile(args.genome + s):
            tblastn.run_makeblastdb(args.genome)
    tblastn_out = tblastn.run_tblastn(args.queries, args.genome)
    tblastn_hits = tblastn.parse_tblastn(tblastn_out)
    unique_hits = tblastn.unique_hits(tblastn_hits)
    ranges = tblastn.join_hits(unique_hits)
    print(ranges)

    # start pulling subsequences
    blocks = subseq.write_subseq(args.genome, ranges)
    # blocks contains tuples of (query, subject, filename)
    for i in blocks:
        exonerate.run_exonerate(i[2])

