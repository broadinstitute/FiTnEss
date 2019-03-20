library(dplyr)
require(graphics)
require(pastecs)
require(ggplot2)

#1. read in raw tally files

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
