#!/usr/bin/env python

## This script finds and reports the position of all TA motifs in a genome
## input:  python TAfinder.py argv[1] > output_file.txt
## argv[1] is the fasta sequence of the organism downloaded form NCBI
## The output file has 2 columns:  1) misc_feature; 2) genome position of the TA site on the Forward (+) strand

#updated to python3, 7/22/25

import sys

fasta=[] #use an array for effecient concatenation of the lines

for line in open(sys.argv[1]):
    # ignore the first line if it starts with >
    if line[0]!=">":
        #now only take the portion of the line up until a tab character
        split = line.split('\t')
        read = split[0].rstrip()
        fasta.append(read)

# join the lines together into a single string, and uppercase it
genome="".join(fasta).upper()

index = genome.find('TA')
while (index != -1):
    print("misc_feature\t%d" %(index + 1))
    index = genome.find('TA', index+2)
