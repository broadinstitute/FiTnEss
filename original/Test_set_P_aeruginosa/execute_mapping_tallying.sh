#building Bowtie genome
bowtie-build ./input/UCBPP-PA14.fasta PA14

#begin running pipeline per sample
echo PA14_M9_rep1
#mapping to genome
bowtie -v 0 -q -m 1 -S PA14 PA14_M9_rep1.fastq PA14_M9_rep1.sam
#tallying
python ./scripts/tally.py ./output/TA_sites.txt ./input/PA14_M9_rep1.sam ./input/genelist_PA14.txt > ./output/PA14_M9_rep1_tally.txt
