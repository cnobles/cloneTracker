#samples need to have PCR replicates removed to reduce the amount of work needed to be done
remove_repeats <- function(sites, mcols.to.keep = NULL){
  sites <- sort(sites)
  
  if(!is.null(mcols.to.keep)){
    keep.pos <- match(mcols.to.keep, names(mcols(sites)))
    cols.dfr <- distinct(data.frame(mcols(sites)[keep.pos]))
    if(nrow(cols.dfr) > 1){
      stop("Cannot have mcols with unique length greater than 1.")
  }}
  
  mcols(sites) <- NULL
  sites.df <- as.data.frame(sites, row.names = NULL)
  sites.df <- distinct(sites.df)
  ranges <- IRanges(start = sites.df$start, end = sites.df$end)
  sites.deamp <- GRanges(seqnames = sites.df$seqnames,
                         ranges = ranges,
                         strand = sites.df$strand,
                         seqinfo = seqinfo(sites))
  if(!is.null(mcols.to.keep)){mcols(sites.deamp) <- cols.dfr}
  sites.deamp
}
