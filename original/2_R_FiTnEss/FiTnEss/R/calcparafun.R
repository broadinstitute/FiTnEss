#calcparafun
#2. calculate parameters

calcparafun<-function(strain,usable_tally_list,save_location,rep_time){

  b=rep_time

  if (strain %in% c("UCBPP","PA14","pa14","ucbpp")){
    strain<-"UCBPP"
    st<-"UCBPP-PA14"
  }else{
    st<-strain
  }

  ####################### Calculating parameters #######################

  print("Starting Phase II: calculating parameters")

  if(st %in% c("PA14","UCBPP","UCBPP-PA14","pa14","ucbpp")){
    st=st
  }else{
    st=st
  }

  print(paste("Now running ",st,sep=""))

  #1. prepare data

  filtered_list<-lapply(usable_tally_list,function(x){
    x=x %>% dplyr::filter(homo==FALSE,non_permissive==FALSE,coreTA==TRUE)
  })

  #2. compare replicates

  repgene_list<-lapply(filtered_list,function(x){
    x=x %>% group_by(Locus.CIA) %>% summarise(Nta=n())
  })
  consistent_genes<-Reduce(function(x, y) inner_join(x, y, by = c("Locus.CIA","Nta")), repgene_list) #5708 genes
  filtered_list<-lapply(filtered_list,function(x){
    x=x %>% dplyr::filter(Locus.CIA %in% unique(consistent_genes$Locus.CIA))
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
    totdata<-oridata %>% dplyr::select(Locus.CIA,n) %>% group_by(Locus.CIA) %>% summarise(gtot=sum(n),Nta=n())
    totdata<-totdata %>% group_by(Nta) %>% mutate(ngene=n()) %>% ungroup()
    nsample<-oridata %>% group_by(Locus.CIA) %>% sample_n(1) %>% dplyr::select(Locus.CIA,n)
    totdata<-totdata %>% left_join(nsample,by="Locus.CIA")
    colnames(totdata)[5]<-"TAsample"
    nsample<-oridata %>% group_by(Locus.CIA) %>% sample_n(1) %>% dplyr::select(Locus.CIA,n)
    tot<-totdata %>% left_join(nsample,by="Locus.CIA")
    a=median(tot$Nta)

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
