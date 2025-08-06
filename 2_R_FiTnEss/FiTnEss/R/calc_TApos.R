
#3. calculate number of TA sites

calc_TApos <- function(tally, genelist) {
  x <- base::merge(tally, genelist, by = c("Chr","Locus.CIA","strain"), all.x = TRUE)
  # x <- tally %>% full_join(genelist,by = c("Locus.CIA","strain"))
  x <- dplyr::filter(x, type == 'CDS') #5858 genes
  x <- dplyr::group_by(x, Chr, Locus.CIA) %>%
    mutate(TA.gene.pos = TA_start - gene_start + 1) %>% #1 would be the TA is at the start of the gene
    mutate(gene.size.CIA = gene_stop - gene_start + 1) %>%
    ungroup()
  x <- as.data.frame(x)
  x$TA.gene.pos[x$strand == '-'] <- x$gene.size.CIA[x$strand == '-']-x$TA.gene.pos[x$strand == '-']+1
  x <- dplyr::mutate(x, TA.gene.percent = TA.gene.pos/gene.size.CIA)
  #x <- x[!duplicated(x$TA_start),]
  x = dplyr::distinct(x) #8/6/25
  x$TAindex <- as.numeric(x$TAindex)
  return(x)
}
