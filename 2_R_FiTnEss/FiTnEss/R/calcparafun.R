#calcparafun
#2. calculate parameters

calcparafun<-function(strain,usable_tally_list,save_location,rep_time){

  b=rep_time

  # if (strain %in% c("UCBPP","PA14","pa14","ucbpp")){
  #   strain<-"UCBPP"
  #   st<-"UCBPP-PA14"
  # }else{
  #   st<-strain
  # }

  ####################### Calculating parameters #######################

  print("Starting Phase II: calculating parameters")

  # if(st %in% c("PA14","UCBPP","UCBPP-PA14","pa14","ucbpp")){
  #   st=st
  # }else{
  #   st=st
  # }

  print(paste("Now running ", strain, sep=""))

  #1. prepare data

  filtered_list<-lapply(usable_tally_list,function(x){
    x=x %>% dplyr::filter(homo==FALSE, non_permissive==FALSE, coreTA==TRUE, remove_plasmid==FALSE)
  })

  #2. compare replicates

  repgene_list<-lapply(filtered_list,function(x){
    x=x %>% group_by(Chr, Locus.CIA) %>% summarise(Nta=n(), .groups = "drop") %>% ungroup()
  })

  #iterative inner join to find genes present in all dataframes in list
  consistent_genes<-Reduce(function(x, y) inner_join(x, y, by = c("Chr","Locus.CIA","Nta")), repgene_list) #5708 genes

  filtered_list<-lapply(filtered_list,function(x){
    #x=x %>% dplyr::filter(Locus.CIA %in% unique(consistent_genes$Locus.CIA))
    #Locus should be unique, but to make it more robust with chromosome info
    x = dplyr::semi_join(x, consistent_genes, by = c("Chr", "Locus.CIA"))
  })

  saveRDS(filtered_list,
          paste(dirname(save_location),
                gsub(".xlsx","_tally.rds",basename(save_location)),sep="/"))
  print("Usable tally files saved")

  #3. calculate parameters

  parameter_list<-lapply(filtered_list,function(x){

    #prepare each replicate
    data=x
    strain<-unique(data$strain)
    oridata<-data
    #Tabulate total number of reads, gtot, and total number of TA sites per gene
    totdata<-oridata %>% dplyr::select(Chr, Locus.CIA, n) %>% group_by(Chr, Locus.CIA) %>% summarise(gtot=sum(n),Nta=n(), .groups = "drop") %>% ungroup()
    #Tabulate how many genes have each # of TA sites
    totdata<-totdata %>% group_by(Nta) %>% mutate(ngene=n()) %>% ungroup()
    #For each gene, pick a random site
    nsample<-oridata %>% group_by(Chr, Locus.CIA) %>% sample_n(1) %>% dplyr::select(Chr, Locus.CIA,n) %>% ungroup()
    totdata<-totdata %>% left_join(nsample,by=c("Chr", "Locus.CIA"))
    colnames(totdata)[6]<-"TAsample"
    #sample again
    nsample<-oridata %>% group_by(Chr, Locus.CIA) %>% sample_n(1) %>% dplyr::select(Chr, Locus.CIA,n) %>% ungroup()
    tot<-totdata %>% left_join(nsample,by=c("Chr", "Locus.CIA"))
    a=median(tot$Nta) #median number of TA sites per gene across the whole genome

    #calculate parameters
    ntalist<-rep(a,b)
    print(ntalist)
    ntareslist<-lapply(ntalist,function(y){
      nta=as.numeric(y)
      samdata=dplyr::filter(tot,Nta==nta)$gtot #376
      y=unlist(optim(c(0.9,log(mean(dplyr::filter(tot,Nta==nta)$gtot)/nta),1),function(pars){my_cvm(samdata,"jointfun",pars[1],pars[2],pars[3],nta)$statistic}))
      print(y)
    })

    ntares<-as.data.frame(do.call("rbind",ntareslist))
    colnames(ntares)<-c("lambda","lp","sigma","cvm","counts.function",
                        "counts.gradient","convergence")
    ntares$Nta<-rep(a,b)

    x=ntares
  })

  names(parameter_list)<-paste("rep",seq_len(length(parameter_list)),sep="_")

  saveRDS(parameter_list,
          paste(dirname(save_location),
                gsub(".xlsx","_parameters.rds",
                     basename(save_location)),sep="/"))
  print("Parameters saved")
  return(parameter_list)

}
