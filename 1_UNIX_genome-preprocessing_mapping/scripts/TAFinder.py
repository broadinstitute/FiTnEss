#!/usr/bin/env python

## This script finds and reports the position of all TA motifs in a genome
## input:  python TAfinder.py argv[1] > output_file.txt
## argv[1] is the fasta sequence of the organism downloaded form NCBI
## The output file has 2 columns:  1) misc_feature; 2) genome position of the TA site on the Forward (+) strand

#updated to python3 & had it save contig names, 7/22/25

import sys

fasta = []
contig_name = None

def process_contig(contig_name, fasta):
    if contig_name and fasta:
        genome = "".join(fasta).upper()
        index = genome.find('TA')
        while index != -1:
            print("%s\t%d" % (contig_name, index + 1))
            index = genome.find('TA', index + 2)

for line in open(sys.argv[1]):
    if line[0] == ">":
        # Process previous contig before starting a new one
        process_contig(contig_name, fasta)
        contig_name = line[1:].strip().split()[0]
        fasta = []
    else:
        split = line.split('\t')
        read = split[0].rstrip()
        fasta.append(read)

# Process the last contig
process_contig(contig_name, fasta)