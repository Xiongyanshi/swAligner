import argparse
import numpy as np
from lib.seqio import Fasta

def theta(a, b):
    """return score of two bases"""
    res = -1
    if a == b:
        res = 1
    return res

class Align:
    def __init__(self, seq1, seq2):
        self.seq1 = '-' + seq1
        self.seq2 = '-' + seq2

    def buildmat(self):
        self.scoremat = np.empty((len(self.seq1), len(self.seq2)), dtype=int)
        self.tracemat = np.empty((len(self.seq1), len(self.seq2)), dtype=int)
        for i, p in enumerate(self.seq1):
            for j, q in enumerate(self.seq2):
                if i == 0:
                    self.scoremat[i, j] == 0
                    self.tracemat[i, j] == 1
                    continue
                if j == 0:
                    self.scoremat[i, j] == 0
                    self.tracemat[i, j] == 2
                    continue
                d = self.scoremat[i-1, j-1] + theta(p, q)
                l = self.scoremat[i, j-1] + theta('-',q)
                u = self.scoremat[i-1, j] + theta(p,'-')
                self.scoremat[i, j] = max(d,l,u,0)
                self.tracemat[i, j] = [d, l, u, 0].index(max(d,l,u,0))


def main():
    descriptionStr = 'Smith-Waterman pairwise alignment(toy) by Yanshi Xiong'
    parser = argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-i','--input', required=True,
                        help="""file name of input fasta,
                              only first 2 records will be read in.""")
    parser.add_argument('-o','--output',
                        help='file name of alignment result(in fasta format)')

    args = parser.parse_args()
    args = vars(args)
    infile = args['input']
    infasta = Fasta(infile)
    outfile = args['output']
    outfasta = Fasta(outfile)

    name1 = infasta.name1
    seq1  = infasta.seq1
    name2 = infasta.name2
    seq2  = infasta.seq2

    align = Align(seq1, seq2)
    align.buildmat()
    print(align.scoremat)
    print(align.tracemat)


if __name__ == '__main__':
    main()
