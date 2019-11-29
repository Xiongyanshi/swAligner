swAligner
> Smith-Waterman Alignment program.

# Usage
## Getting start
```
# quick help
python swAligner.py -h

# run demo
python swAligner.py

# input sequences from file
python swAligner.py -i demo/input1.fa

# input sequences from copy and paste
python swAligner.py -r atactttcg,acccttttcgaa

# write aligned reads into file
python swAligner.py -i demo/input2.fa -o input2.aligned.fa

# use another base substitution matrix
python swAligner.py -i demo/input2.fa -m lib/dna_sub1.mat

# use another gap penalty
python swAligner.py -i demo/input2.fa -g lib/dna_gap.2.txt

```
## Advanced usage
### Prepare your own substitution matrix
substitution matrix is a tab delimited table, first character of first row doesn't matter, only int will be accepted (default, lib/dna_sub.default.mat):
```
\	a	t	c	g
a	3	-3	-3	-3
t	-3	3	-3	-3
c	-3	-3	3	-3
g	-3	-3	-3	3
```
For instance, if you prefer more similarity between purines or pyrimidines, you may wanna use this one (lib/dna_sub.1.mat) :
```
\	a	t	c	g
a	3	-3	-3	1
t	-3	3	1	-3
c	-3	1	3	-3
g	1	-3	-3	3
```
### Another gap penalty scheme
Default gap penalty configuration file is "lib/dna_gap.default.txt" :
```
open=-2
extend=-2
```
It means opening a gap will -2, and extending a gap -2 as well. In other words, there is no difference between opening a gap and extending a gap. As such, "-" match to "A" means -2, "---" match to "ACG" means -6

But, you may prefer opening a gap cost much more than extending a gap, this will make alignment program trying to find continuous body, for example in lib/dna_gap.1.txt:
```
open=-5
extend=-2
```
In this way, "-" match to "A" means -5, "---" match to "ACG" means -9 (not -15).

For instance:
```
$ python swAligner.py  -i demo/input2.fa -g lib/dna_gap.default.txt

TACGGGCCCGCTAC------
|| |  | | |||
TA-G--C-C-CTATCGGTCA

$ python swAligner.py  -i demo/input2.fa -g lib/dna_gap.1.txt

TACGGGCCCGCTAC------
     || | |||
---TAGC-C-CTATCGGTCA

```

# Algorithm
## How it works
The Smith-Waterman algorithm performs local sequence alignment, a derived version of Needleman-Wunsch algorithm. Based on match/mismatch/gap penalty scheme, a m*n matrix is generated and followed by tracing back from the highest score to the zero cell. Two matrix is built in this python program: "scoremat", a matrix of the scores calculated, and "tracemat", a matrix of which way the value in scoremat is caculated. The alignment solution is generated by the tracing back the direction in "tracemat", and represented in "pathcode". The match/mismatch score (substitution matrix) and gap penally are both configurable.

## Module required
```
re, numpy, argparse
```


