# to execute this file in terminal, type: sh execute_TAsite_preprocessing.sh
# modify fasta filenames on lines 13, 16, 19, and 20. Example shown is for "UCBPP-PA14".fasta

# required files in folder:
# 1. TAFinder.py
# 2. nonpermissive_TAFinder.py
# 3. 50bp_TA_sequence_extractor.py
# 4. homology_finder.py
# 5. your genome fasta file
# 6. your genome fasta file with sequence in reverse complement. Many free tools exist such as: http://reverse-complement.com

#TA site finder. Finds all TA sites in genome with fasta file input
python ./scripts/TAFinder.py ./input/UCBPP-PA14.fasta > ./output/TA_sites.txt

#nonpermissive TA site finder. Finds all TA sites that are non-permissive to mariner transposon insertion in genome with fasta file input
python ./scripts/nonpermissive_TAFinder.py ./input/UCBPP-PA14.fasta > ./output/nonpermissive_TA_sites.txt

#homology finder. Finds all TA sites that occur more than once in a genome with 2 fasta files as inputs
python ./scripts/50bp_TA_sequence_extractor.py ./input/UCBPP-PA14.fasta > 50bp_TA_sequences.txt
python ./scripts/50bp_TA_sequence_extractor.py ./input/UCBPP-PA14_reverse_complement.fasta > 50bp_TA_sequences_reverse.txt
tail -r 50bp_TA_sequences_reverse.txt > temp.txt
cut -f1 temp.txt > temp2.txt
cut -f2 50bp_TA_sequences.txt | paste temp2.txt - > temp3.txt
cat temp3.txt 50bp_TA_sequences.txt > 50bp_TA_sequences_combined.txt
python ./scripts/homology_finder.py 50bp_TA_sequences_combined.txt > ./output/homologous_TA_sites.txt

#removal of temp files
rm temp.txt
rm temp2.txt
rm temp3.txt
rm 50bp_TA_sequences.txt
rm 50bp_TA_sequences_reverse.txt
rm 50bp_TA_sequences_combined.txt