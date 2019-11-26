import argparse
from lib.seqio import Seq

def main():
    parser = argparse.ArgumentParser(description='Smith-Waterman pairwise alignment (toy) by Yanshi Xiong')
    parser.add_argument('-i','--input', required=True,
                        help='file name of input fasta for alignment, only first 2 sequences will be read in.')
    parser.add_argument('-o','--output', required=True,
                        help='file name of alignment result in fasta format')

    args = parser.parse_args()
    inputfile = args['input']
    outfile = args['output']

    name1, seq1, name2, seq2 = Seq(inputfile)
    print(args)
    print(name1)
    print(seq1)
    print(name2)
    print(seq2)

if __name__ == '__main__':
    main()
