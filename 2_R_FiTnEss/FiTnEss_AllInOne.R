#FiTnEss_AllInOne
#Oct15_2018_TS_pipeline_function
library(dplyr)
library(tidyr)
library(fBasics)
library(goftest)
library(openxlsx)
# library(pastecs)
library(scales)
library(stats)


FiTnEss_Run<-function(strain,
                      supporting_functions,
                      file_location,
                      permissive_file,
                      homologous_file,
                      gene_file,
                      save_location,
                      repeat_time=3){ #pre-defined as using Nta=10, and replicate for 5 times for each replicate

  source(supporting_functions)

  ## 1. Usable tally file preparation

  usable_tally_list<-self_tallyprepfun(strain,file_location,permissive_file,homologous_file,gene_file)

  ## 2. Calculating parameters

  parameter_list<-calcparafun(strain,usable_tally_list,save_location,rep_time=repeat_time)

  ## 3. Call Essentials

  result_list<-callessfun(file_location,usable_tally_list,parameter_list)

  ## 4. save final results

  write.xlsx(result_list, file = save_location)
  print("Final results saved, finished running.")

}












































