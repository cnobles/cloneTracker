db_to_granges <- function(dfr_from_db, keep.additional.columns = FALSE){
  dfr <- dfr_from_db
  ranges <- IRanges(start = ifelse(dfr$strand == "+", dfr$position, dfr$breakpoint),
                    end = ifelse(dfr$strand == "+", dfr$breakpoint, dfr$position))
  gr <- GRanges(seqnames = dfr$chr,
                ranges = ranges,
                strand = dfr$strand)
  
  GTSP.list <- strsplit(dfr$sampleName, split="-")
  mcols <- data.frame("sampleName" = dfr$sampleName,
                      "GTSP" = sapply(1:length(gr), function(i){GTSP.list[[i]][1]})
                      stringAsFactors = FALSE)
  
  if(keep.additional.columns){
    std.columns <- c("sampleName", "position", "chr", "strand", "breakpoint")
    are.there <- match(std.columns, colnames(dfr))
    add.cols <- grep(TRUE, is.na(match(names(dfr), names(dfr[,are.there]))))
    cols <- data.frame(dfr[, add.cols])
    colnames(cols) <- colnames(dfr[add.cols])
    mcols <- cbind(mcols, cols)
  }
  
  mcols(gr) <- mcols
  gr
}
