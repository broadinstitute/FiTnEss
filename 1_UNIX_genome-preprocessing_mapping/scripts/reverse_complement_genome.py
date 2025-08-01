#!/usr/bin/env python3

#Create reverse complement of a genome fasta file using Biopython

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

fna_file_path = sys.argv[1]
output_file_path = sys.argv[2]

def process_fasta_with_biopython(input_file, output_file):
    """Process FASTA file and write reverse complements using Biopython."""
    try:
        sequences_processed = 0

        with open(output_file, 'w') as output_handle:
            # Parse input FASTA file
            for record in SeqIO.parse(input_file, "fasta"):
                # Create reverse complement using Biopython
                rev_comp_seq = record.seq.reverse_complement()

                # Create new record with modified description
                new_id = f"{record.id}"
                new_description = f"{record.description}"

                # Create new SeqRecord
                from Bio.SeqRecord import SeqRecord
                rev_comp_record = SeqRecord(
                    rev_comp_seq,
                    id=new_id,
                    description=new_description
                )

                # Write to output file
                SeqIO.write(rev_comp_record, output_handle, "fasta")
                sequences_processed += 1

        return sequences_processed

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None
    except Exception as e:
        print(f"Error processing file: {e}")
        return None


def validate_sequences(input_file):
    """Validate and show information about input sequences."""
    try:
        sequences = list(SeqIO.parse(input_file, "fasta"))

        if not sequences:
            print("No sequences found in the input file.")
            return False

        print(f"Found {len(sequences)} sequence(s):")

        for i, record in enumerate(sequences, 1):
            seq_type = "DNA" if set(str(record.seq).upper()) <= set("ATGCNRYSWKMBDHV") else "Unknown"
            print(f"  {i}. {record.id} - Length: {len(record.seq)} bp - Type: {seq_type}")

            # Show first 50 bases as preview
            preview = str(record.seq)[:50]
            if len(record.seq) > 50:
                preview += "..."
            print(f"     Preview: {preview}")

        return True

    except Exception as e:
        print(f"Error validating sequences: {e}")
        return False
    

#Run functions to create reverse complement .fna
if(not os.path.exists(output_file_path)):
  # Validate input sequences
  if validate_sequences(fna_file_path):
    print("-" * 50)
    print("Processing sequences...")

    # Process the file
    sequences_processed = process_fasta_with_biopython(fna_file_path, output_file_path)

    if sequences_processed is not None:
        print(f"Successfully processed {sequences_processed} sequence(s)")
        print(f"Reverse complement sequences written to: {output_file_path}")
        print("Done!")
    else:
        print("Processing failed.")
  else:
    print("Validation failed - check your input file.")
else:
   print(output_file_path + "already exists. Not overwritting.")
