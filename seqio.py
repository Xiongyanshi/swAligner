import re

class Seq:
    def __init__(self, filename):
        self.__filecontent = open(filename, 'r').read()
        try:
            name1, seq1, name2, seq2 = re.match(r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)\n>.+', self.__filecontent).groups()
        except AttributeError:
            name1, seq1, name2, seq2 = re.match(r'^>(.+?)\n([\w\W]+?)\n>(.+?)\n([\w\W]+)', self.__filecontent).groups()
        seq1 = ''.join(seq1.splitlines())
        seq2 = ''.join(seq2.splitlines())

        self.name1 = name1
        self.seq1  = seq1
        self.name2 = name2
        self.seq2  = seq2

        self.read1 = '>%s\n%s' % (name1, seq1)
        self.read2 = '>%s\n%s' % (name2, seq2)

    def write(self, outfilename):
        with open(outfilename, 'w') as handle:
            handle.write('>%s\n%s\n>%s\n%s' % (self.name1, self.seq1, self.name2, self.seq2))

