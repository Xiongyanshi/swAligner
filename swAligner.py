import argparse
from seqio import Seq

def main():
    parser = argparse.ArgumentParser(description='Smith-Waterman pairwise alignment (toy) by Yanshi Xiong')
    parser.add_argument('-i','--input', required=True, help='file name of input fasta for alignment, only first 2 sequences will be read in.')
    parser.add_argument('-o','--output', required=True, help='file name of alignment result in fasta format')

    args = parser.parse_args()
    print(args)

if __name__ == '__main__':
    main()
