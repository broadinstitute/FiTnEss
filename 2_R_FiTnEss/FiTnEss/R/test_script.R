#Test each function with the new input files
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("fBasics", quietly = TRUE)) {
  install.packages("fBasics")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}


library(dplyr)
library(openxlsx)
library(fBasics)
library(tidyr)

strain = "PA14"
working_dir = "/Volumes/idi_hunglabusers/jbagnall/klebs_tnseq/test/Test_set_P_aeruginosa/"
setwd(working_dir)
file_location = "./sample_data/PA14_M9_rep1_test_tally.txt"
permissive_file = "./TAsite_info/nonpermissive_TA_sites.txt"
homologous_file = "./TAsite_info/homologous_TA_sites.txt"
gff_file = "./genome_info/GCA_000014625.1.gff"
save_location = "./test_results_250805.xlsx"
gff_name_tag = "locus_tag"
repeat_time = 3

FiTnEss_Run(strain = strain,
                      file_location = file_location,
                      permissive_file = permissive_file,
                      homologous_file = homologous_file,
                      gene_file = gff_file,
                      save_location = save_location,
                      gff_name_tag = gff_name_tag, #tag to define gene name in gff file
                      repeat_time = repeat_time)

# usable_tally_list = self_tallyprepfun(strain = strain, file_location = file_location, permissive_file = permissive_file, homologous_file = homologous_file, gene_file = gff_file, gff_name_tag = gff_name_tag)

#fit function
# parameter_list<-calcparafun(strain,usable_tally_list,save_location,rep_time=repeat_time)

# result_list<-callessfun(file_location,usable_tally_list,parameter_list)
