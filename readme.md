swAligner
> Smith-Waterman Alignment program.

# Usage
## getting start
```
# quick help
python swAligner.py -h

# run demo
python swAligner.py

# input sequence from file
python swAligner.py -i demo/input1.fa

# input sequence from raw paste
python swAligner.py -r atactttcg,acccttttcgaa

# write aligned reads into file
python swAligner.py -i demo/input2.fa -o input2.aligned.fa

# use another base substitution matrix
python swAligner.py -i demo/input2.fa -m lib/dna_sub1.mat
```
## Advanced usage:
### prepare your own substitution matrix
substitution matrix is a tab delimited table, first character of first row doesn't matter, only int will be passed (the default one in lib/dna_default.mat):
```
\	a	t	c	g
a	3	-3	-3	-3
t	-3	3	-3	-3
c	-3	-3	3	-3
g	-3	-3	-3	3
```
For instance, if you prefer more similarity between purines or pyrimidines, you may wanna use this one( in lib/dna_sub1.mat) :
```
\	a	t	c	g
a	3	-3	-3	1
t	-3	3	1	-3
c	-3	1	3	-3
g	1	-3	-3	3
```

# Algorithm
## how it works

## module required
```
re, numpy, argparse
```

## Smith-Waterman alignment
