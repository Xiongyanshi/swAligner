swAligner
> Smith-Waterman Alignment program.

# Usage
## getting start
```
# run demo
python swAligner.py -i demo/input1.fa

# quick help
python swAligner.py -h

# write aligne reads into file
python swAligner.py -i demo/input2.fa -o input2.aligned.fa

# use another base substitution matrix
python swAligner.py -i demo/input1.fa -m lib/dna_default.mat
```
## import it as a module

# Algorithm
## module required
re,numpy, argparse
## Smith-Waterman alignment
