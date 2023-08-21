#! /usr/bin/python3


# courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
def parse_fasta(path):
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


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
        self.queries = [self.query]  # for merging, if multi queries


class genome_range():
    """
    simple referenceable data class for genomic spans
    """
    def __init__(self, input=[]):
        self.query = input[0]
        self.subject = input[1]
        self.start = input[2]
        self.end = input[3]
        self.queries = [self.query]
        self.file = ""
