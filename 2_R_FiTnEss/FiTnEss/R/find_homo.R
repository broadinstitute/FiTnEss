
#2. denotate TAs that are homologous
#8/5/25 adjusted to allow for multiple contigs

find_homo<-function(tally,homofile){
  homofile2<-homofile %>%
    group_by(V2,V3) %>% #include contig
    summarise(side=n(), .groups = "drop") %>%
    ungroup() %>%
    select(V2,V3,side)
  colnames(homofile2)<-c("Chr", "TA_start","side")
  tally = dplyr::mutate(tally, Chr_TAstart = paste(Chr, TA_start, sep = ":"))
  homofile2 = dplyr::mutate(homofile2, Chr_TAstart = paste(Chr, TA_start, sep = ":"))
  tally$homo<-FALSE
  tally$homo[which(tally$Chr_TAstart %in% homofile2$Chr_TAstart)]<-TRUE
  tally<-tally %>% left_join(homofile2,by=c("Chr", "TA_start", "Chr_TAstart"))
  tally$side[which(is.na(tally$side))]<-0
  return(tally)
}
