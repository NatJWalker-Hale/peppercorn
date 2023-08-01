#! /usr/bin/python3


import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def write_subseq(dbf, ranges: list[tuple]) -> list[tuple]:
    coord_dict = {}
    for x in ranges:
        try:
            coord_dict[x[1]].append((x[0], x[2], x[3]))
        except KeyError:
            coord_dict[x[1]] = [(x[0], x[2], x[3])]
    # print(coord_dict)
    out = []
    for rec in SeqIO.parse(dbf, "fasta"):
        if rec.id in coord_dict.keys():
            for coords in coord_dict[rec.id]:  # tuples
                output = "_".join([dbf, rec.id, str(coords[1]),
                                   str(coords[2])])
                newid = rec.id + "_" + str(coords[1]) + "_" + str(coords[2])
                # print((coords[1]-1, coords[2]-1))
                newrec = SeqRecord(
                    Seq(rec.seq[coords[1]-1:coords[2]-1]),
                    id=newid,
                    description="",
                    )
                SeqIO.write(newrec, output, "fasta")
                out.append((coords[0], rec.id, output))
    return out


# rewrite this for lists and break into two functions
def redundant_ranges(coord_dict: dict, buffer=300) -> dict:
    """
    if ranges are within buffer, combine
    two cases
    ----------------
        -----------------

    overlap is r1 end - r2 start

        -----------------
    ----------------

    overlap is r2 end - r1 start

    if overlap > buffer, merge
    """
    out = {}
    for k, v in coord_dict.items():
        for x, y in coord_dict.items():
            if k == x:  # same contig
                overlap = False
                # overlap = 0
                for r1 in v:
                    for r2 in y:
                        if r1[0] < r2[0] and r1[1] > r2[0]:
                            overlap_len = r1[1] - r2[0]
                            if overlap_len > buffer:
                                overlap = True
                                try:
                                    out[k].append((r1[0], r2[1]))
                                except KeyError:
                                    out[k] = [(r1[0], r2[1])]
                            y.remove(r2)
                        elif r2[0] < r1[0] and r2[1] > r1[0]:
                            overlap_len = r2[1] - r1[0]
                            if overlap_len > buffer:
                                overlap = True
                                try:
                                    out[k].append((r2[0], r1[1]))
                                except KeyError:
                                    out[k] = [(r2[0], r1[1])]
                            y.remove(r2)
                        else:  # no overlap
                            continue
                        # overlap = min(r1[1], r2[1]) - max(r1[0], r2[0])
                        # if overlap > 0:
                        #     newbound = (min(r1[0], r2[0]), max(r1[1], r2[1]))
                        # try:
                        #     out[k].append(newbound)
                        # except KeyError:
                        #     out[k] = [newbound]
                        # v.remove(r1)
                # if overlap == 0:
                if not overlap:
                    try:
                        out[k].append(r1)
                    except KeyError:
                        out[k] = [r1]
    return out


def range_to_coord_dict(ranges: list[tuple]) -> dict:
    """
    takes output like pad_range and converts to coord_dict as expected by
    subseq utilities
    """
    out = {}
    for r in ranges:
        try:
            out[r[1]].append((r[2], r[3]))
        except KeyError:
            out[r[1]] = [(r[2], r[3])]
    return out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("database", help="FASTA-formatted database to \
                        extract subsequences")
    parser.add_argument("id", help="which sequence to extract from")
    parser.add_argument("range", help="space-separated start and end range \
                        to extract, 1-indexed, e.g. 24252 24353",
                        type=int, nargs=2)
    args = parser.parse_args()

    coord_dict = {"OW204023.1": [(24000, 25000), (24500, 25500)],
                  "OW204024.1": [(24000, 25000), (24500, 25500)]}

    coord_dict = redundant_ranges(coord_dict, 300)
    print(coord_dict)

    files = write_subseq(args.database, coord_dict)
    print(files)
