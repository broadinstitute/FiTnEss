#Script to align reads and tally reads per TA site per gene

###################Input variables--modify these######################################################################
#test manuscript data
data_path='~/test/Test_set_P_aeruginosa_v2/' #location of genome_info and TAsite_info folders
script_path='~/FiTnEss_JB/1_UNIX_genome-preprocessing_mapping/scripts/' #location of scripts
temp_path='~/temp/250807/' #where to output temp files
fna_file='GCA_000014625.1_ASM1462v1_genomic.fna'
gff_file='GCA_000014625.1.gff' #GFF file with gene annotations, used to create gene list file
fastq_file='M91_PA14.fastq.gz' #fastq file with reads to map
tally_output_name="PA14_M9_rep1_tally.txt" #output file name for tally results, usually has strain_media

#############No need to modify below unless changing number of allowed mismatches######################################
#building Bowtie genome
base_fna_name=$(basename "${fna_file}" .fna) #base name of fasta file, used to create reverse complement fasta file name
index_file="${base_fna_name}_index"
bowtie-build "${data_path}genome_info/${fna_file}" "${data_path}genome_info/${index_file}"

#begin running pipeline per sample
#mapping to genome, exact match
base_fastq_name=$(basename "${fastq_file}" .fastq.gz)
bowtie -v 0 -q -m 1 -S "${data_path}genome_info/${index_file}" "${data_path}sample_data/${fastq_file}" "${data_path}sample_data/${base_fastq_name}.sam"
#Allowing 1 mismatch
#bowtie -v 1 -q -m 1 -S "${data_path}genome_info/${index_file}" "${data_path}sample_data/${fastq_file}" "${data_path}sample_data/${base_fastq_name}.sam"

#tallying"${script_path}reverse_complement_genome.py"
base_gff_name=$(basename "${gff_file}" .gff) #base name of GFF file, used to create gene list file name
gene_list_name="${base_gff_name}_gene_list.txt"
python "${script_path}tally.py" "${data_path}TAsite_info/TA_sites.txt" "${data_path}sample_data/${base_fastq_name}.sam" "${data_path}genome_info/${gene_list_name}" "${data_path}genome_info/${fna_file}" > "${data_path}sample_data/${tally_output_name}"
