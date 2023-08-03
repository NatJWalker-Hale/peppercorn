#! /usr/bin/python3


import os
import sys
import utils
import argparse
import tblastn
import subseq
import exonerate
import augustus


def split_query(queryf, queries=None):
    queryseqs = dict([x for x in utils.parse_fasta(queryf)])
    if queries is not None:
        if len(queries) > 1:
            with open("multiquery.pep.fa", "w") as outf:
                for k, v in queryseqs.items():
                    if k in queries:
                        outf.write(">%s\n" % k)
                        outf.write("%s\n" % v)
            return "multiquery.pep.fa"
        else:
            with open(queries[0] + ".pep.fa", "w") as outf:
                for k, v in queryseqs.items():
                    if k in queries:
                        outf.write(">%s\n" % k)
                        outf.write("%s\n" % v)
            return queries[0] + ".pep.fa"
    else:
        with open(k + ".pep.fa", "w") as outf:
            for k, v in queryseqs.items():
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
    parser.add_argument("-i", "--intronlength", help="Expected maximum intron \
                        length (if unsure, longer better)", type=int,
                        default=3000)
    parser.add_argument("-t", "--threads", help="Number of threads", type=int,
                        default=1)
    args = parser.parse_args()

    # run tblastn
    blastdbsuf = [".ndb", ".nhr", ".nin", ".nog", ".nos", ".not",
                  ".nsq", ".ntf", ".nto"]
    for s in blastdbsuf:
        if not os.path.isfile(args.genome + s):
            tblastn.run_makeblastdb(args.genome)
    tblastn_out = tblastn.run_tblastn(args.queries, args.genome, args.threads)
    tblastn_hits = tblastn.parse_tblastn(tblastn_out)
    hit_lens = [x.sendorder - x.sstartorder for x in tblastn_hits]
    mean_hit_len = round(sum(hit_lens) / len(hit_lens))
    # unique_hits = tblastn.unique_hits(tblastn_hits)
    ranges = tblastn.hits_to_range(tblastn.join_hits(tblastn_hits))
    print([(x.query, x.subject,
            x.start, x.end,
            x.file) for x in ranges])
    tblastn.pad_range(ranges, length=args.intronlength + mean_hit_len)
    # start pulling subsequences
    subseq.write_subseq(args.genome, ranges)
    print([(x.query, x.subject,
            x.start, x.end,
            x.file) for x in ranges])
    for r in ranges:
        print(r.query)
        print(r.queries)
        qf = split_query(args.queries, r.queries)
        out, hintsf, res = exonerate.run_exonerate(r.file, qf)
        if not res:
            sys.stderr.write("No exonerate result for %s, discarding\n"
                             % r.file)
            os.remove(r.file)
            os.remove(hintsf)
            # os.remove(res)
            ranges.remove(r)
        else:
            augout = augustus.run_augustus(r.file, hintsf)
            print(augustus.parse_augustus(augout))

        # stuff with augustus
