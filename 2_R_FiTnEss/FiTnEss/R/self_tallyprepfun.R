#self-prepare supporting files
#self_tallyprepfun
#tallyprepfun


#Function to extract attribute from gff
extract_tag <- function(attr_string, tag) {
  # Create case-insensitive regex
  pattern <- paste0("(?i)", tag, "=[^;]+")  # (?i) turns on case-insensitive matching

  match <- regmatches(attr_string, regexpr(pattern, attr_string, perl = TRUE))

  if (length(match) == 0 || match == "") {
    return(NA)  # return NA if tag not found
  }

  sub(paste0("(?i)", tag, "="), "", match, perl = TRUE)
}


#1. usable tally file preparation

self_tallyprepfun<-function(strain, file_location, permissive_file, homologous_file, gene_file, gff_name_tag = "locus_tag", remove_multicopy_plasmid_names){

  #genefile<-read.delim(gene_file, header = FALSE, fill = TRUE, sep = '\t') #original script, not parsing gff correctly
  genefile <- read.table(gene_file, sep = "\t", header = FALSE, quote = "", comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)

  ####################### Usable tally file preparation #######################

  print("Starting Phase I: preparing usable tally files")

  usable_tally_list<-list()
  file_location_list<-unlist(strsplit(file_location,","))

  #1. load in supporting files

  # if (strain %in% c("UCBPP","PA14","pa14","ucbpp")){
  #   strain<-"UCBPP" #why are these different?
  #   st<-"PA14"
  # }else{
  #   st<-strain
  # }
  st = strain #why is there a distinction
  nm<-paste(st,"_support",sep="")

  nonPermissiveTA<-read.delim(permissive_file, header = FALSE, fill = TRUE)
  colnames(nonPermissiveTA) = c("Seq", "Chr", "TA_start")
  nonPermissiveTA = dplyr::mutate(nonPermissiveTA, Chr_TAstart = paste(Chr, TA_start, sep = ":"))

  homo<-read.delim(homologous_file, header = FALSE, fill = TRUE)

  genelist<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V1,V9,V4,V5)
  genelist$V9 <- sapply(genelist$V9, extract_tag, tag = gff_name_tag)

  geneinfo<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V1,V3,V7,V9)
  geneinfo$V9 <- sapply(geneinfo$V9, extract_tag, tag = gff_name_tag)

  cluster2<-genefile %>% dplyr::filter(V3=="CDS") %>% dplyr::select(V1,V9)
  colnames(cluster2) = c("Chr", "V9")
  cluster2$V92<-cluster2$V9
  cluster2$Locus.CIA<-sapply(cluster2$V9, extract_tag, tag = gff_name_tag)
  cluster2$desc<-sapply(cluster2$V9, extract_tag, tag = "gene") #could also do product

  cluster2$strain<-strain
  cluster2<-cluster2 %>% dplyr::select(Chr, Locus.CIA,strain,desc)

  colnames(genelist)<-c("Chr","Locus.CIA","gene_start","gene_stop")
  # genelist$strain<-st
  genelist$strain = strain
  colnames(geneinfo)<-c("Chr", "type","strand","Locus.CIA")
  genelist<-genelist %>% dplyr::full_join(geneinfo,by=c("Chr", "Locus.CIA")) %>% dplyr::distinct()
  genelist<-genelist %>% dplyr::left_join(cluster2,by=c("Chr", "Locus.CIA","strain"))
  genelist<-genelist %>% dplyr::select(Chr,Locus.CIA,strain,type,gene_start,gene_stop,strand,desc) %>% dplyr::distinct()

  # if (strain=="UCBPP"){
  #   genelist$strain<-"UCBPP"
  # }else{
  #   genelist$strain<-genelist$strain
  # }
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
    # unique_map_tally2$strain<-strain

    #b) Import raw tally files and annotate sites with multialignment:
    unique_map_tally <- calc_TApos(unique_map_tally2, genelist)
    ## removed many TA sites
    unique_map_tally <- dplyr::filter(unique_map_tally, type == 'CDS')
    tally.w <- unique_map_tally

    #c) Find sites affected by non-permissive bias:

    # tally.w$non_permissive <- (tally.w$TA_start %in% nonPermissiveTA$V2)
    tally.w$non_permissive <- (tally.w$Chr_TAstart %in% nonPermissiveTA$Chr_TAstart)
    table(tally.w$non_permissive) #TRUE is number of non-permissive TA sites

    #d) Denote TAs for edge trimming

    tally.w <- denote_coreTA(tally.w, 50)
    ## denote TAs as TRUE it is not in the first and last 50bp of the gene
    table(tally.w$coreTA) #false is number of TA sites found near edges

    #mark any in a plasmid or contig that should be removed
    #initialize to false
    tally.w$remove_plasmid = FALSE
    if(all(!is.na(remove_multicopy_plasmid_names))){
      tally.w = dplyr::mutate(tally.w, remove_plasmid = tolower(Chr) %in% tolower(remove_multicopy_plasmid_names))
      if(all(tally.w$remove_plasmid==FALSE)){
        print(paste0("Warning: no genes found in contigs ", paste(remove_multicopy_plasmid_names, collapse = ",")))
      }
    }

    #e) process to list for downstream analysis
    x=tally.w

  })
  return(usable_tally_list)
}
