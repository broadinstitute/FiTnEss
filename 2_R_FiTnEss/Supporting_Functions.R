#Supporting_Functions
library(dplyr)
library(tidyr)
library(goftest)
library(fBasics)
library(openxlsx)

######################
#  Data preparation  #
######################

import_raw_tally <- function(filename, CIA = TRUE) {
  df <- read.delim(filename, header = FALSE) #100654
  if (CIA) {
    df <- df[-grep('^IG',df$V4),]
  } else {
    df <- df
  }
  colnames(df) <- c('Chr','TA_start','TA_end','Locus.CIA','plus_n','minus_n')
  df <- mutate(df, n = plus_n + minus_n)
  df$strain <- strsplit(filename,'_')[[1]][2]
  ## ---IMPORTANT! pay attention to the sequence of file names
  df$Media <- strsplit(filename,'_')[[1]][1]
  #df$Media <- substr(med[i], 1, nchar(med[i])-1)
  rownames(df) <- NULL
  df$TAindex <- rownames(df)
  return(df)
}

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

calc_TApos <- function(tally, genelist) {
  x <- base::merge(tally, genelist, by = c('Locus.CIA',"strain"), all.x = TRUE)
  # x <- tally %>% full_join(genelist,by = c("Locus.CIA","strain"))
  x <- dplyr::filter(x, type == 'CDS') #5858 genes
  x <- dplyr::group_by(x, Locus.CIA) %>%
    mutate(TA.gene.pos = TA_start - gene_start + 1) %>%
    mutate(gene.size.CIA = gene_stop - gene_start + 1)
  x <- as.data.frame(x)
  x$TA.gene.pos[x$strand == '-'] <- x$gene.size.CIA[x$strand == '-']-x$TA.gene.pos[x$strand == '-']+1
  x <- mutate(x, TA.gene.percent = TA.gene.pos/gene.size.CIA)
  x <- x[!duplicated(x$TA_start),]
  x$TAindex <- as.numeric(x$TAindex)
  return(x)
}

denote_coreTA <- function(tally, bp = 50) {
  tally$coreTA <- ((tally$TA.gene.pos > bp) &
                     (tally$TA.gene.pos < (tally$gene.size.CIA-bp + 1)))
  return(tally)
}

##############################
#  Distribution preparation  #
##############################

my_cvm = function (x, null = "punif", ..., nullname){
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if (is.character(null))
    nulltext <- null
  if (missing(nullname) || is.null(nullname)) {
    reco <- goftest::recogniseCdf(nulltext)
    nullname <- if (!is.null(reco))
      reco
    else paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if (is.function(null))
    null
  else if (is.character(null))
    get(null, mode = "function")
  else stop("Argument 'null' should be a function, or the name of a function")
  U <- F0(x, ...)

  if (any(U < 0 | U > 1)){
    omega2<-Inf
    out <- list(statistic = omega2, p.value = NA,
                method = "Coerce into Inf when null function outside [0,1]",
                data.name = xname)
    class(out) <- "htest"
    return(out)
  }
  U <- sort(U)
  x <- sort(x)
  xmn = mean(x)
  xsd = sd(x)
  k <- seq_len(n)

  omega2 <- 1/(12 * n) + sum(fBasics::Heaviside(k,0.25*n)*(U - (2 * k - 1)/(2 * n))^2)
  ## -- changed to 0.5
  PVAL <- goftest::pCvM(omega2, n = n, lower.tail = FALSE)
  names(omega2) <- "omega2"
  METHOD <- c("Cramer-von Mises test of goodness-of-fit",
              paste("Null hypothesis:",nullname))
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if (length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for (i in seq_along(parnames))
      pard[i] <- paste(parnames[i],"=", paste(pars[[i]], collapse = " "))
    pard <- paste("with", ngettext(length(pard), "parameter","parameters"),
                  "  ", paste(pard, collapse = ", "))
    METHOD <- c(METHOD, pard)
  }
  out <- list(statistic = omega2, p.value = PVAL, method = METHOD,
              data.name = xname)
  class(out) <- "htest"
  return(out)
}

nblgnfun2<-function(q,lp,sigma,nta){
  lbd1<-1/(rlnorm(100000,meanlog=lp,sdlog=sigma))
  ndist<-rnbinom(100000,nta,lbd1)
  pvec=sapply(q,function(x){sum(x>=ndist,na.rm = TRUE)/length(ndist)})
  return(pvec)
}

lownbfun3<-function(q,nta){
  pvec=pnbinom(q,nta,0.7)
  return(pvec)
}

jointfun<-function(q,lambda,lp,sigma,nta){
  pvec=lambda*nblgnfun2(q,lp,sigma,nta)+(1-lambda)*lownbfun3(q,nta)
  return(pvec)
}

########################
#  Prepare Tally File  #
########################

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

##########################
#  Calculate Parameters  #
##########################

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
    st="UCBPP"
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

##########################
#  Call essential genes  #
##########################

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
    totdata<-oridata %>% dplyr::select(Locus.CIA,n) %>%
      group_by(Locus.CIA) %>% summarise(gtot=sum(n),Nta=n())
    totdata<-totdata %>% group_by(Nta) %>% mutate(ngene=n()) %>% ungroup()
    
    para=parameter_list
    lp=para[[x]]$lp[which(para[[x]]$cvm==min(para[[x]]$cvm))]
    sigma=para[[x]]$sigma[which(para[[x]]$cvm==min(para[[x]]$cvm))]
    
    ngene<-totdata %>% ungroup() %>% group_by(Nta) %>% summarise(ngene=n())
    
    plist<-lapply(seq_len(length(ngene$Nta)),function(x){
      lbd1<-1/rlnorm(10000, meanlog = lp, sdlog = sigma)
      Zg1<-rnbinom(1000000,ngene$Nta[x],lbd1)
      Zg2<-rnbinom(1000000,ngene$Nta[x],0.7)
      subtot<-totdata %>% dplyr::filter(Nta==ngene$Nta[x])
      subtot2<-subtot %>% rowwise() %>%
        mutate(pv1=length(which(Zg1<=gtot))/length(Zg1),
               pv2=length(which(Zg2>=gtot))/length(Zg2)) #pv2: the larger the better
      subtot3<-subtot2
      x=subtot3
    })
    names(plist)<-paste("Nta",ngene$Nta,sep="_")
    pdata1<-do.call("rbind",plist)
    
    #adjusted p-value
    
    pdata1<-pdata1 %>% dplyr::select(Locus.CIA:pv1)
    colnames(pdata1)[5]<-"pvalue"
    #FWER
    pdata1$padj<-p.adjust(pdata1$pvalue)
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















