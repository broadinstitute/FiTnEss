#!/usr/bin/env python3

## This script uses the upstream or downstream extracted sequences and looks for replicate sequences within the file
## input example:  python homology_finder.py 50bp_downstream_TA_sequences.txt > downstream_50bp_homologous_TAsites.txt
## The output file has 2 columns:  1) the replicated surrounding TA sequence; 2) genome position of the TA site on the Forward (+) strand

#python3 compatible, 7/22/25

import sys

homology={}

for line in open(sys.argv[1]):
    split = line.split('\t')
    sequence = split[0]
    list = homology.get(sequence)
    if list == None:
        list = []
        homology[sequence] = list
    list.append(line)

for sequence, list in homology.items():
    if len(list) > 1:
        for line in list:
            print(line.rstrip())

