#!/usr/bin/env python3

## This script finds and reports the position of all TA motifs in a genome and the downstream sequence
## input:  python downstream_TA_sequence_extractor.py filename.fasta > downstream_TA_sequences.txt
## The output file has 2 columns:  1) TA + 48 nucleotides downstream of TA; 2) genome position of the TA site on the Forward (+) strand
#python3 compatible, adjusted to keep track of multiple contigs 7/22/25

#!/usr/bin/env python3

## This script finds and reports the position of all TA motifs in a genome and the downstream sequence
## input:  python downstream_TA_sequence_extractor.py filename.fasta > downstream_TA_sequences.txt
## The output file has 3 columns:  
## 1) TA + 48 nucleotides downstream of TA
## 2) contig name
## 3) position of the TA site within the contig (1-based)

import sys

def process_contig(sequence, contig_name):
    """Find all TA sites in a contig and print their sequences and positions."""
    index = sequence.find('TA')
    while index != -1:
        ##if index + 50 <= len(sequence):  # don't require 50bp available since then the revrse compelment may not have the same number of sites
        downstream_sequence = sequence[index:index+50]
        print(f"{downstream_sequence}\t{contig_name}\t{index + 1}")
        index = sequence.find('TA', index+2)

def read_fasta(filename):
    """Read FASTA file and process each contig separately."""
    current_contig = None
    current_seq = []
    
    for line in open(filename):
        if line[0] == ">":
            # Process previous contig if it exists
            if current_contig and current_seq:
                sequence = "".join(current_seq).upper()
                process_contig(sequence, current_contig)
                current_seq = []
            current_contig = line[1:].strip().split()[0]
        else:
            split = line.split('\t')
            read = split[0].rstrip()
            current_seq.append(read)
    
    # Process the last contig
    if current_contig and current_seq:
        sequence = "".join(current_seq).upper()
        process_contig(sequence, current_contig)

if __name__ == "__main__":
    read_fasta(sys.argv[1])
