library(dplyr)
require(pastecs)

#2. denotate TAs that are homologous

find_homo<-function(tally,homofile){
  homofile2<-homofile %>% group_by(V2) %>%
    summarise(side=n()) %>% dplyr::select(V2,side)
  colnames(homofile2)<-c("TA_start","side")
  tally$homo<-FALSE
  tally$homo[which(tally$TA_start %in% homofile2$TA_start)]<-TRUE
  tally<-tally %>% left_join(homofile2,by="TA_start")
  tally$side[which(is.na(tally$side))]<-0
  return(tally)
}
