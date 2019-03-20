#!/usr/bin/env python

## This script finds and reports the position of all TA motifs in a genome and the downstream sequence
## input:  python downstream_TA_sequence_extractor.py filename.fasta > downstream_TA_sequences.txt
## The output file has 2 columns:  1) TA + 48 nucleotides downstream of TA; 2) genome position of the TA site on the Forward (+) strand

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
    downstream_sequence = genome[index:index+50]
    print "{}\t{}".format(downstream_sequence,(index + 1))
    index = genome.find('TA', index+2)
