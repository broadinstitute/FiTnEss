library(dplyr)
require(graphics)
require(pastecs)
require(ggplot2)

#3. calculate number of TA sites

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
