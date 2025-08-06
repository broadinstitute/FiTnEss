#callessfun
#3. calculate essential genes

callessfun<-function(file_location,usable_tally_list,parameter_list){

  file_location_list<-unlist(strsplit(file_location,","))
  filtered_list<-lapply(usable_tally_list,function(x){
    x=x %>% dplyr::filter(homo==FALSE,non_permissive==FALSE,coreTA==TRUE)
  })

  ####################### Call Essentials #######################

  print("Starting Phase III: calling Essentials")

  result_list<-list()
  repvec<-seq_len(length(unique(file_location_list)))

  result_list<-lapply(seq_len(length(unique(file_location_list))),function(x){
    data=filtered_list[[x]]
    strain<-unique(data$strain)
    oridata<-data
    totdata<-oridata %>% dplyr::select(Chr, Locus.CIA, n) %>%
      group_by(Chr, Locus.CIA) %>% summarise(gtot=sum(n),Nta=n(), .groups = "drop") %>% ungroup()
    # totdata<-totdata %>% group_by(Nta) %>% mutate(ngene=n()) %>% ungroup()

    para=parameter_list
    lp=para[[x]]$lp[which(para[[x]]$cvm==min(para[[x]]$cvm))]
    sigma=para[[x]]$sigma[which(para[[x]]$cvm==min(para[[x]]$cvm))]

    ngenedf<-totdata %>% ungroup() %>% group_by(Nta) %>% summarise(ngene=n()) %>% ungroup()

    plist<-lapply(seq_len(length(ngenedf$Nta)),function(y){
      lbd1<-1/rlnorm(10000, meanlog = lp, sdlog = sigma)
      Zg1<-rnbinom(1000000,ngenedf$Nta[y],lbd1)
      Zg2<-rnbinom(1000000,ngenedf$Nta[y],0.7)
      subtot<-totdata %>% dplyr::filter(Nta==ngenedf$Nta[y])
      subtot2<-subtot %>% rowwise() %>%
        mutate(pv1=length(which(Zg1<=gtot))/length(Zg1),
               pv2=length(which(Zg2>=gtot))/length(Zg2)) #pv2: the larger the better
      subtot3<-subtot2
      y=subtot3
    })
    names(plist)<-paste("Nta",ngenedf$Nta,sep="_")
    pdata1<-do.call("rbind",plist)

    #adjusted p-value
    pdata1<-pdata1 %>% dplyr::select(Chr:pv1)
    colnames(pdata1)[5]<-"pvalue"
    #FWER
    pdata1$padj<-p.adjust(pdata1$pvalue, "holm") #used default, which should be Holm
    pdata1$Ess_fwer<-"NE_fwer"
    pdata1$Ess_fwer[which(pdata1$padj<0.05)]<-"E_fwer"
    #FDR
    pdata1$pfdr<-p.adjust(pdata1$pvalue,"fdr")
    pdata1$Ess_fdr<-"NE_fdr"
    pdata1$Ess_fdr[which(pdata1$pfdr<0.05)]<-"E_fdr"
    x=pdata1

  })

  print("Finished running, save...")
  names(result_list)<-names(parameter_list)
  return(result_list)
}
