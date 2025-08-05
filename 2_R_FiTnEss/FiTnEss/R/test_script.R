#Test each function with the new input files

strain = "PA14"
working_dir = "/Volumes/idi_hunglabusers/jbagnall/klebs_tnseq/test/Test_set_P_aeruginosa/"
setwd(working_dir)
file_location = "./sample_data/PA14_M9_rep1_test_tally.txt"
permissive_file = "./TAsite_info/nonpermissive_TA_sites.txt"
homologous_file = "./TAsite_info/homologous_TA_sites.txt"
gene_file = "./genome_info/GCA_000014625.1.gff"
save_location = "./test_results_250805.xlsx"
gff_name_tag = "locus_tag"
repeat_time = 3

library(dplyr)
library(tidyr)

self_tallyprepfun(strain = strain, file_location = file_location, permissive_file = permissive_file, homologous_file = homologous_file, gene_file = gff_file, gff_name_tag = gff_name_tag)



temp1 <- read.table(gff_file, sep = "\t", header = FALSE, quote = "", comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
