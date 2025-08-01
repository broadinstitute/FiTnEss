#!/usr/bin/env python3

## This script finds and reports the position of all TA motifs in a genome
## input:  python TAfinder.py argv[1] > output_file.txt
## argv[1] is the fasta sequence of the organism downloaded form NCBI
## The output file has 2 columns:  1) misc_feature; 2) genome position of the TA site on the Forward (+) strand

#updated to python3, 7/22/25

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
    fasta=[] #use an array for effecient concatenation of the lines

    for line in open(filename):
        # ignore the first line if it starts with >
        if line[0] != ">":
            #now only take the portion of the line up until a tab character
            split = line.split('\t')
            read = split[0].rstrip()
            fasta.append(read)

    # join the lines together into a single string, and uppercase it
    return "".join(fasta).upper()


def find_without_overlap(genome):
    ###
    # Use this function if there should not be overlap.
    # As soon as we find a matching sequence, that sequence can no longer be used in another one.
    # The order that we define the patterns to search for is important.
    ###
    pattern = re.compile("|".join(values))
    for match in pattern.finditer(genome):
        print_it(match.group(), match.start())
        # print match.group(), (match.start() + 1)

def find_with_overlap( genome ):
    ###
    # Use this function if there can be overlap in the sequence indexes.
    # For example, given 'ABCDEFG', and searching for 'ABCD' and 'DEFG', the 'D' can be used in both.
    ###
    for value in values:
        # print value
        index = find_next(genome, value, 0)
        while (index != -1):
            # print "%s\t%d" %(value, (index + 1))
            print_it(value, index)
            index = find_next(genome, value, index+2)#genome.find('CGATAACG', index+2)

def find_next( genome, value, previous_index ):
    return genome.find(value, previous_index)

def print_it(value, index):
    print("%s\t%d" % (value, (index + 1 + 3)))

genome = read_data(filename)

# print '------: w/o overlap :------'
# find_without_overlap(genome)

# print '------: w/ overlap :------'
find_with_overlap(genome)
