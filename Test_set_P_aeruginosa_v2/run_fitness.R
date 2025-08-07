#Test each function with the new input files
#Load packages
packages <- c("devtools", "dplyr", "fBasics", "goftest", "openxlsx", "scales", "stats", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if (!require("FiTnEss", character.only = TRUE)) {
  devtools::install_github("broadinstitute/FiTnEss", subdir = "2_R_FiTnEss/FiTnEss")
}

library(FiTnEss)

Packages <- c("dplyr","fBasics","goftest","openxlsx","scales","stats","tidyr")
lapply(Packages, library, character.only = TRUE)

require(FiTnEss)


#Run test file (Brad's data)
strain = "PA14"
working_dir = "~/test/Test_set_P_aeruginosa_v2/" #set your working directory, which has all the files and to which it will save results
setwd(working_dir)
file_location = "./sample_data/PA14_M9_rep1_tally.txt"
permissive_file = "./TAsite_info/nonpermissive_TA_sites.txt"
homologous_file = "./TAsite_info/homologous_TA_sites.txt"
gff_file = "./genome_info/GCA_000014625.1.gff"
save_location = "./test_results_250807.xlsx"
gff_name_tag = "locus_tag"
repeat_time = 5

FiTnEss_Run(strain = strain,
            file_location = file_location,
            permissive_file = permissive_file,
            homologous_file = homologous_file,
            gene_file = gff_file,
            save_location = save_location,
            gff_name_tag = gff_name_tag, #tag to define gene name in gff file
            repeat_time = repeat_time)
