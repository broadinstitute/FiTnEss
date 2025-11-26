#FiTnEss_function
#fitnessRun_function
#Oct23_2018_FiTnEss_separate
#Oct23_2018_FiTnEss

#4. function to calculate parameters using Nta=median(Nta)

FiTnEss_Run<-function(strain,
                      file_location,
                      permissive_file,
                      homologous_file,
                      gene_file,
                      save_location,
                      gff_name_tag, #tag to define gene name in gff file
                      remove_multicopy_plasmid_names = NA, #add plasmid names matching the gff, to remove them
                      repeat_time=3){

  ## 1. Usable tally file preparation

  usable_tally_list<-self_tallyprepfun(strain = strain, file_location = file_location, permissive_file = permissive_file,
                     homologous_file = homologous_file, gene_file = gene_file, gff_name_tag = gff_name_tag,
                     remove_multicopy_plasmid_names = remove_multicopy_plasmid_names)
  ## 2. Calculating parameters

  parameter_list<-calcparafun(strain = strain, usable_tally_list = usable_tally_list, save_location = save_location, rep_time=repeat_time)

  ## 3. Call Essentials

  result_list<-callessfun(file_location = file_location, usable_tally_list = usable_tally_list, parameter_list = parameter_list)

  ## 4. save final results

  openxlsx::write.xlsx(result_list, file = save_location)
  print("Final results saved, finished running.")

}

