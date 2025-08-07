# to execute this file in terminal, type: sh execute_TAsite_preprocessing.sh
# modify input paths for data and scripts

# required files in scripts folder:
# 1. TAFinder.py
# 2. nonpermissive_TAFinder.py
# 3. 50bp_TA_sequence_extractor.py
# 4. homology_finder.py
# 5. reverse_complement_genome.py
# 6. make_gene_list.py

conda activate py3.7.11

#Input variables, 7/22/25
#test from brad's data
data_path='/idi/hunglabusers/jbagnall/klebs_tnseq/test/Test_set_P_aeruginosa_v2/' #location of genome_info and TAsite_info folders
script_path='/home/unix/jbagnall/git/FiTnEss_JB/1_UNIX_genome-preprocessing_mapping/scripts/' #location of scripts
temp_path='/broad/hptmp/jbagnall/temp/250807/' #where to output temp files
fna_file='GCA_000014625.1_ASM1462v1_genomic.fna'
gff_file='GCA_000014625.1.gff' #GFF file with gene annotations, used to create gene list file
gene_name_tag='locus_tag' #tag in GFF file to use for gene names, usually 'locus' or 'locus_tag' or 'ID''

#Make sure necessary subdirectories exist
mkdir -p "${temp_path}"
mkdir -p "${data_path}genome_info" "${data_path}TAsite_info"

#Create reverse complement of genome fasta file if it does not already exist
base_fna_name=$(basename "${fna_file}" .fna) #base name of fasta file, used to create reverse complement fasta file name
fna_reverse_file="${base_fna_name}_reverse_complement.fna"
python "${script_path}reverse_complement_genome.py" "${data_path}genome_info/${fna_file}" "${data_path}genome_info/${fna_reverse_file}"  

#Create gene list file from GFF file if it does not already exist (output filen name is automatically generated)
python "${script_path}make_gene_list.py" "${data_path}genome_info/${gff_file}" "$gene_name_tag"

#TA site finder. Finds all TA sites in genome with fasta file input
python "${script_path}TAFinder.py" "${data_path}genome_info/${fna_file}" > "${data_path}"/TAsite_info/TA_sites.txt

#nonpermissive TA site finder. Finds all TA sites that are non-permissive to mariner transposon insertion in genome with fasta file input
python "${script_path}nonpermissive_TAFinder.py" "${data_path}genome_info/${fna_file}" > "${data_path}TAsite_info/nonpermissive_TA_sites.txt"

#homology finder. Finds all TA sites that occur more than once in a genome with 2 fasta files as inputs
python "${script_path}50bp_TA_sequence_extractor.py" "${data_path}genome_info/${fna_file}" > "${temp_path}50bp_TA_sequences.txt"
python "${script_path}50bp_TA_sequence_extractor.py" "${data_path}genome_info/${fna_reverse_file}" > "${temp_path}50bp_TA_sequences_reverse.txt"
#tail -r 50bp_TA_sequences_reverse.txt > temp.txt , original line threw error


# Split into separate files by contig
awk '{
    print $0 > "'${temp_path}'seq_" $2 ".txt"
}' "${temp_path}50bp_TA_sequences_reverse.txt"

awk '{
    print $0 > "'${temp_path}'pos_" $2 ".txt"
}' "${temp_path}50bp_TA_sequences.txt"

# Reverse and combine matching contigs (only specific columns)
for contig in $(awk '{print $2}' "${temp_path}50bp_TA_sequences.txt" | uniq); do
    tac "${temp_path}seq_${contig}.txt" | cut -f1 > "${temp_path}seq_${contig}_rev.txt"
    cut -f2,3 "${temp_path}pos_${contig}.txt" > "${temp_path}pos_${contig}_cut.txt"
    paste "${temp_path}seq_${contig}_rev.txt" "${temp_path}pos_${contig}_cut.txt" > "${temp_path}combined_${contig}.txt"
done

# Concatenate all contig files
cat "${temp_path}combined_"*.txt > "${temp_path}temp3.txt"
cat "${temp_path}temp3.txt" "${temp_path}50bp_TA_sequences.txt" > "${temp_path}50bp_TA_sequences_combined.txt"

# Find homologous TA sites
python "${script_path}homology_finder.py" "${temp_path}50bp_TA_sequences_combined.txt" > "${data_path}TAsite_info/homologous_TA_sites.txt"

# Clean up temporary files
rm "${temp_path}"seq_*.txt "${temp_path}"pos_*.txt "${temp_path}"combined_*.txt "${temp_path}temp3.txt"
