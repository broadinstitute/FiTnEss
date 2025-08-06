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
                      repeat_time=3){ #pre-defined as using Nta=10, and replicate for 5 times for each replicate

  ## 1. Usable tally file preparation

  usable_tally_list<-self_tallyprepfun(strain,file_location,permissive_file,homologous_file,gene_file)

  ## 2. Calculating parameters

  parameter_list<-calcparafun(strain,usable_tally_list,save_location,rep_time=repeat_time)

  ## 3. Call Essentials

  result_list<-callessfun(file_location,usable_tally_list,parameter_list)

  ## 4. save final results

  openxlsx::write.xlsx(result_list, file = save_location)
  print("Final results saved, finished running.")

}

