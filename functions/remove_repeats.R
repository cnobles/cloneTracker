#samples need to have PCR replicates removed to reduce the amount of work needed to be done
remove_repeats <- function(sites){
  sites.df <- as.data.frame(sites)
  sites.df <- distinct(sites.df)
  ranges <- IRanges(start = sites.df$start, end = sites.df$end)
  sites.deamp <- GRanges(seqnames = sites.df$seqnames,
                         ranges = ranges,
                         strand = sites.df$strand,
                         seqinfo = seqinfo(sites))
  sites.deamp
}
