#!/usr/bin/env python3

## This script finds and reports the position of all TA motifs in a genome
## input:  python TAfinder.py argv[1] > output_file.txt
## argv[1] is the fasta sequence of the organism downloaded form NCBI
## The output file has 2 columns:  1) misc_feature; 2) genome position of the TA site on the Forward (+) strand

#updated to python3, and to handle multiple contigs 7/22/25

import sys
import re

filename = sys.argv[1]
values = [
    'CGATAACG',
    'CGTTAACG',
    'CGGTAACG',
    'CGCTAACG',
    'CGATATCG',
    'CGTTATCG',
    'CGGTATCG',
    'CGCTATCG',
    'CGATAGCG',
    'CGTTAGCG',
    'CGGTAGCG',
    'CGCTAGCG',
    'CGATACCG',
    'CGTTACCG',
    'CGGTACCG',
    'CGCTACCG',
    'CGATAACC',
    'CGTTAACC',
    'CGGTAACC',
    'CGCTAACC',
    'CGATATCC',
    'CGTTATCC',
    'CGGTATCC',
    'CGCTATCC',
    'CGATAGCC',
    'CGTTAGCC',
    'CGGTAGCC',
    'CGCTAGCC',
    'CGATACCC',
    'CGTTACCC',
    'CGGTACCC',
    'CGCTACCC',
    'GGATAACG',
    'GGTTAACG',
    'GGGTAACG',
    'GGCTAACG',
    'GGATATCG',
    'GGTTATCG',
    'GGGTATCG',
    'GGCTATCG',
    'GGATAGCG',
    'GGTTAGCG',
    'GGGTAGCG',
    'GGCTAGCG',
    'GGATACCG',
    'GGTTACCG',
    'GGGTACCG',
    'GGCTACCG',
    'GGATAACC',
    'GGTTAACC',
    'GGGTAACC',
    'GGCTAACC',
    'GGATATCC',
    'GGTTATCC',
    'GGGTATCC',
    'GGCTATCC',
    'GGATAGCC',
    'GGTTAGCC',
    'GGGTAGCC',
    'GGCTAGCC',
    'GGATACCC',
    'GGTTACCC',
    'GGGTACCC',
    'GGCTACCC'
]


def read_data(filename):
    """Read FASTA file and return dictionary of contig sequences."""
    contigs = {}
    current_contig = None
    current_seq = []

    for line in open(filename):
        if line[0] == ">":
            if current_contig:
                contigs[current_contig] = "".join(current_seq).upper()
                current_seq = []
            current_contig = line[1:].strip().split()[0]
        else:
            split = line.split('\t')
            read = split[0].rstrip()
            current_seq.append(read)
    
    if current_contig:  # Don't forget the last contig
        contigs[current_contig] = "".join(current_seq).upper()
    
    return contigs

def find_with_overlap(contigs):
    """Find motifs in each contig separately."""
    for contig_name, sequence in contigs.items():
        for value in values:
            index = find_next(sequence, value, 0)
            while (index != -1):
                print("%s\t%s\t%d" % (value, contig_name, (index + 1 + 3)))
                index = find_next(sequence, value, index+2)

def find_next(sequence, value, previous_index):
    return sequence.find(value, previous_index)

# Remove or comment out the old print_it function as it's no longer needed
# def print_it(value, index):
#     print("%s\t%d" % (value, (index + 1 + 3)))

contigs = read_data(filename)
find_with_overlap(contigs)