
#1. read in raw tally files

import_raw_tally <- function(filename, CIA = TRUE) {
  df <- read.delim(filename, header = FALSE) #100654
  if (CIA) {
    df <- df[-grep('^IG',df$V4),]
  } else {
    df <- df
  }
  colnames(df) <- c('Chr','TA_start','TA_end','Locus.CIA','plus_n','minus_n')
  df <- dplyr::mutate(df, n = plus_n + minus_n)

  #Do we really want to assume the file ame is strain_media?
  base_filename = basename(filename)
  df$strain <- strsplit(base_filename,'_')[[1]][1]
  ## ---IMPORTANT! pay attention to the sequence of file names
  df$Media <- strsplit(base_filename,'_')[[1]][2]
  #df$Media <- substr(med[i], 1, nchar(med[i])-1)
  rownames(df) <- NULL
  df$TAindex <- rownames(df)
  return(df)
}
