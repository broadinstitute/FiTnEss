#building Bowtie genome
bowtie-build ./genome_info/your_genome.fasta ./genome_info/your_genome

#begin running pipeline per sample
echo your_sample
#mapping to genome
bowtie -v 0 -q -m 1 -S ./genome_info/your_genome ./sample_data/your_sample.fastq ./sample_data/your_sample.sam
#tallying
python ./scripts/tally.py ./TAsite_info/TA_sites.txt ./sample_data/your_sample.sam ./genome_info/your_genelist.txt > ./sample_data/your_sample_tally.txt
