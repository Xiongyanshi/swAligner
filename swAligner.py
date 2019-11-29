#! /usr/bin/env python

"""
A Smith-Waterman alignment program.
by Yanshi Xiong.
"""

import argparse
import numpy as np
from lib.seqio import Fasta
from lib.readmat import readmat
from lib.readgap import readgap


def tailn(l, x):
    """
    return the number of continuous element x from the right side of l
    tail([0,1,2,3,4,3,3,3,3], 3) = 4
    """
    count = 0
    for i in list(l)[::-1]:
        if i != x:
            break
        count += 1
    return count

class Align:
    def __init__(self, seq1, seq2, submat, gap_open, gap_extend):
        self.seq1 = '-' + seq1.upper()
        self.seq2 = '-' + seq2.upper()
        self.submat = submat
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        Align.__buildmat(self)
        Align.__traceback(self)
        Align.__makeprintable(self)

    def __buildmat(self):
        self.scoremat = np.empty([len(self.seq1), len(self.seq2)],
                                 dtype=np.int32)
        self.tracemat = np.empty([len(self.seq1), len(self.seq2)],
                                 dtype=np.int8)
        for i, p in enumerate(self.seq1):
            for j, q in enumerate(self.seq2):
                if i == 0:
                    self.scoremat[i, j] = 0
                    self.tracemat[i, j] = 1
                    continue
                if j == 0:
                    self.scoremat[i, j] = 0
                    self.tracemat[i, j] = 2
                    continue

                # from diagnoal
                d = self.scoremat[i-1, j-1] + self.submat[p][q]

                # from left side
                lgaplen = tailn(self.tracemat[i,:j], 1)
                lgap_penalty = lgaplen * self.gap_open + self.gap_extend
                l = self.scoremat[i, j-1] + lgap_penalty + self.submat[p][q]

                # from up side
                ugaplen = tailn(self.tracemat[:i,j], 2)
                ugap_penalty = ugaplen * self.gap_open + self.gap_extend
                u = self.scoremat[i-1, j] + ugap_penalty + self.submat[p][q]

                self.scoremat[i, j] = max(d,l,u,0)
                self.tracemat[i, j] = [d, l, u, 0].index(max(d,l,u,0))

    def __traceback(self):
        self.maxscore = self.scoremat.max()
        self.maxscorei = np.unravel_index(self.scoremat.argmax(),
                                          self.scoremat.shape)
        self.pathcode = ''
        i, j = self.maxscorei   # this is where we start to traceback
        while self.scoremat[i,j] != 0:
            #print('[%s,%s]>' % (i, j), end='')
            direction = str(self.tracemat[i,j])
            self.pathcode = direction + self.pathcode
            if direction == '0':
                i -= 1
                j -= 1
            if direction == '1':
                j -= 1
            if direction == '2':
                i -= 1
            if direction == '3':   # this is not in plan.
                print('wrong')
                print(self.scoremat)
                print(self.tracemat)
                print(i, j)
        self.start = [i,j]
        self.end = self.maxscorei

    def __makeprintable(self):
        pathcode = self.pathcode
        i, j = self.start  # start of aligned base
        i = i+1
        j = j+1

        # the unaligned seq at left side
        leftsize = max(i, j) - 1                   # -1 : cut '-' place
        top    = '{:>{leftsize}}'.format(self.seq1[1:i], leftsize=leftsize)
        middle = ' ' * leftsize
        bottom = '{:>{leftsize}}'.format(self.seq2[1:j], leftsize=leftsize)

        # connect with aligned body
        for path in pathcode:
            if path == '0':
                top    += self.seq1[i]
                middle += '|' if self.seq1[i] == self.seq2[j] else '*'
                bottom += self.seq2[j]
                i += 1
                j += 1
            if path == '1':
                top    += '-'
                middle += ' '
                bottom += self.seq2[j]
                j += 1
            if path == '2':
                top    += self.seq1[i]
                middle += ' '
                bottom += '-'
                i += 1

        # tail with unaligned sequence at right side
        m, n = self.end
        m += 1
        n += 1
        rightsize = max(len(self.seq1) - m, len(self.seq2) - n)
        top    += '{:<{rightsize}}'.format(self.seq1[m:], rightsize=rightsize)
        middle += ' ' * rightsize
        bottom += '{:<{rightsize}}'.format(self.seq2[n:], rightsize=rightsize)

        top = top.replace(' ', '-')
        bottom = bottom.replace(' ', '-')
        res = "\n%s\n%s\n%s\n" % (top, middle, bottom)
        self.alignedseq1 = top
        self.alignedseq2 = bottom
        self.printable = res

    def print_align(self):
        print(self.printable)

    def print_scoremat(self):
        print(self.scoremat)

    def print_tracemat(self):
        print(self.tracemat)

    def __repr__(self):
        return self.printable


def main():
    descriptionStr = 'Smith-Waterman pairwise alignment(toy) by Yanshi Xiong'
    parser = argparse.ArgumentParser(description=descriptionStr)
    parser.add_argument('-i','--input', required=False,
                        default='./demo/input1.fa',
                        help='file name of input fasta, parse first 2 records')
    parser.add_argument('-o','--output', required=False,
                        default=None,
                        help='file name of alignment result (fasta)')
    parser.add_argument('-r', '--reads', required=False,
                        default=None,
                        help='two reads seperated by "," (skip open file)')
    parser.add_argument('-m','--submat', required=False,
                        default='./lib/dna_sub.default.mat',
                        help='file name of base substitution matrix')
    parser.add_argument('-g', '--gap', required=False,
                        default='./lib/dna_gap.default.txt',
                        help='gap penalty plan')

    args = vars(parser.parse_args())

    # substitution matrix and gap penalty are must.
    submat     = readmat(args['submat'])
    gap_open   = readgap(args['gap'])['gap_open']
    gap_extend = readgap(args['gap'])['gap_extend']

    # get input
    if args['reads']:
        seq1, seq2 = args['reads'].split(',')
    else:
        infasta = Fasta(args['input'])
        seq1, seq2 = infasta.seq1, infasta.seq2

    # ready to go
    align = Align(seq1, seq2, submat, gap_open, gap_extend)
    align.print_align()
    #align.print_scoremat()
    #align.print_tracemat()

    # if writing to file
    if args['output']:
        outfile = args['output']
        outfasta = Fasta()
        outfasta.from_string('>%s\n%s\n>%s\n%s\n' % (
                               infasta.name1 + '_aligned',
                               align.alignedseq1,
                               infasta.name2 + '_aligned',
                               align.alignedseq2))
        outfasta.to_file(outfile)

    return 0

if __name__ == '__main__':
    main()
