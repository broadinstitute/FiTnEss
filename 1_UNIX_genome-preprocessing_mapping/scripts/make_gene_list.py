#!/usr/bin/env python3

#Create gene list file using Biopython
#input: GFF file with gene annotations, which tag to use for gene name (e.g. locus_tag), and output file path
#output: gene list file with columns: chromosome/contig, gene_id, blank, blank, blank, start, end

import re
import sys
import os
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    from Bio import SeqIO
    from Bio.Seq import Seq

gff_file_path = sys.argv[1]
tag_name = sys.argv[2]


def create_gene_list_from_gff(gff_file, output_file=None):
    #Create a gene list file from a GFF file.
    #Args:
    #    gff_file (str): Path to input GFF file
    #    output_file (str): Path to output gene list file (optional)
    #Returns:
    #    str: Path to created gene list file
    
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(gff_file))[0]
        output_file = f"{base_name}_genelist.txt"

    genes_found = 0

    try:
        with open(gff_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue

                # Split GFF line (tab-delimited)
                fields = line.split('\t')

                # GFF format: seqname, source, feature, start, end, score, strand, frame, attribute
                if len(fields) < 9:
                    continue

                seqname = fields[0]     # chromosome/contig
                feature = fields[2]     # feature type
                start = fields[3]       # start position
                end = fields[4]         # end position
                attributes = fields[8]  # attributes field

                # Only process gene features
                if feature.lower() != 'cds':
                    continue

                # Extract gene ID from attributes
                gene_id = extract_gene_id(attributes, tag=tag_name)

                if gene_id:
                    # Format: chr, gene_id, blank, blank, blank, start, end
                    outfile.write(f"{seqname}\t{gene_id}\t\t\t\t{start}\t{end}\n")
                    genes_found += 1

        print(f"Created gene list: {output_file}")
        print(f"Found {genes_found} genes")
        return output_file

    except FileNotFoundError:
        print(f"Error: GFF file '{gff_file}' not found.")
        return None
    except Exception as e:
        print(f"Error processing GFF file: {e}")
        return None
    

def extract_gene_id(attributes, tag="locus_tag"):
    """
    Extract a tag value (default: locus_tag) from a GFF attributes string.
    
    Parameters:
        attributes (str): The semicolon-separated attribute field from a GFF.
        tag (str): The attribute tag to extract (e.g., 'locus_tag', 'ID').
        
    Returns:
        str or None: The value associated with the tag, or None if not found.
    """
    # Make case-insensitive search pattern
    pattern = re.compile(rf'{re.escape(tag)}=([^;]+)', re.IGNORECASE)
    match = pattern.search(attributes)
    
    if match:
        return match.group(1).strip('"')

    return None


#Run function to make gene list
#Generate output filename automatically
base_name = os.path.splitext(gff_file_path)[0]
output_gff = f"{base_name}_gene_list.txt"

# Check if the file exists and is empty before skipping
if os.path.exists(output_gff) and os.path.getsize(output_gff) > 0:
  print(output_gff + " already exists and is not empty. Not overwritting.")
else:
  create_gene_list_from_gff(gff_file_path, output_gff)