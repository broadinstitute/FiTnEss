#self-prepare supporting files
#self_tallyprepfun
#tallyprepfun

#1. usable tally file preparation

self_tallyprepfun<-function(strain,file_location,permissive_file,homologous_file,gene_file){
  
  genefile<-read.delim(gene_file, header = FALSE, fill = TRUE)
  
  ####################### Usable tally file preparation #######################
  
  print("Starting Phase I: preparing usable tally files")
  
  usable_tally_list<-list()
  file_location_list<-unlist(strsplit(file_location,","))
  
  #1. load in supporting files
  
  if (strain %in% c("UCBPP","PA14","pa14","ucbpp")){
    strain<-"UCBPP"
    st<-"PA14"
  }else{
    st<-strain
  }
  nm<-paste(st,"_support",sep="")
  
  nonPermissiveTA<-read.delim(permissive_file, header = FALSE, fill = TRUE)
  homo<-read.delim(homologous_file, header = FALSE, fill = TRUE)
  genelist<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V1,V9,V4,V5)
  genelist$V9<-sapply(genelist$V9,function(x){
    a=as.character(x)
    x=gsub("locus=","",strsplit(a,";")[[1]][3])
  })
  geneinfo<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V3,V7,V9)
  geneinfo$V9<-sapply(geneinfo$V9,function(x){
    a=as.character(x)
    x=gsub("locus=","",strsplit(a,";")[[1]][3])
  })
  
  cluster2<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V9)
  cluster2$V92<-cluster2$V9
  cluster2$Locus.CIA<-sapply(cluster2$V9,function(x){
    a=as.character(x)
    x=gsub("locus=","",strsplit(a,";")[[1]][3])
  })
  cluster2$desc<-sapply(cluster2$V9,function(x){
    a=as.character(x)
    x=gsub("name=","",strsplit(a,";")[[1]][4])
  })
  cluster2$strain<-strain
  cluster2<-cluster2 %>% dplyr::select(Locus.CIA,strain,desc)
  
  colnames(genelist)<-c("V1","Locus.CIA","gene_start","gene_stop")
  genelist$strain<-st
  colnames(geneinfo)<-c("type","strand","Locus.CIA")
  genelist<-genelist %>% full_join(geneinfo,by="Locus.CIA")
  genelist<-genelist %>% left_join(cluster2,by=c("Locus.CIA","strain"))
  genelist<-genelist %>% dplyr::select(Locus.CIA,strain,type,gene_start,gene_stop,strand,desc)
  
  if (strain=="UCBPP"){
    genelist$strain<-"UCBPP"
  }else{
    genelist$strain<-genelist$strain
  }
  print("Finished loading supporting files")
  
  #2. prepare usable tally files
  
  usable_tally_list<-lapply(file_location_list,function(x){
    unique_file=x
    print(unique_file)
    
    #a) Import raw tally files and annotate sites with homologous:
    
    unique_map_tally <- import_raw_tally(unique_file,CIA=FALSE)
    unique_map_tally<-find_homo(unique_map_tally,homo)
    unique_map_tally2<-unique_map_tally
    unique_map_tally2$Locus.CIA<-gsub("IG_","",unique_map_tally2$Locus.CIA)
    unique_map_tally2$strain<-strain
    
    #b) Import raw tally files and annotate sites with multialignment:
    
    unique_map_tally <- calc_TApos(unique_map_tally2, genelist)
    ## removed many TA sites
    unique_map_tally <- dplyr::filter(unique_map_tally, type == 'CDS')
    tally.w <- unique_map_tally
    
    #c) Find sites affected by non-permissive bias:
    
    tally.w$non_permissive <- (tally.w$TA_start %in% nonPermissiveTA$V2)
    table(tally.w$non_permissive) #TRUE is number of non-permissive TA sites
    
    #d) Denote TAs for edge trimming
    
    tally.w <- denote_coreTA(tally.w, 50)
    ## denote TAs as TRUE it is not in the first and last 50bp of the gene
    table(tally.w$coreTA) #false is number of TA sites found near edges
    
    #e) process to list for downstream analysis
    x=tally.w
    
  })
  return(usable_tally_list)
}
