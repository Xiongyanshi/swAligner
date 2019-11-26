import re

class FastaIO:
    def __init__(self):
        pass

    def parse_String(self, stringContent):
        """parse sequence file, get two reads"""
        try:
            # more than 2 > in fasta file?
            name1, seq1, name2, seq2 = re.match(
                                 r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)\n>.+',
                                 stringContent).groups()
        except AttributeError:
            # no just 2 > records in fasta file, safe to parse until file end
            name1, seq1, name2, seq2 = re.match(
                                 r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)',
                                 stringContent).groups()
        seq1 = ''.join(seq1.splitlines())
        seq2 = ''.join(seq2.splitlines())

        self.name1 = name1
        self.seq1  = seq1
        self.name2 = name2
        self.seq2  = seq2

        self.read1 = '>%s\n%s' % (name1, seq1)
        self.read2 = '>%s\n%s' % (name2, seq2)

        return self

    def parse_File(self, filename):
        return FastaIO.parse_String(open(filename).read())

    def to_file(self, outfilename):
        with open(outfilename, 'w') as handle:
            handle.write('>%s\n%s\n>%s\n%s' % (self.name1, self.seq1,
                                               self.name2, self.seq2))

