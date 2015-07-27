db_to_GRanges <- function(dfr_from_db){
  dfr <- dfr_from_db
  ranges <- IRanges(start = ifelse(dfr$strand == "+", dfr$position, dfr$breakpoint),
                    end = ifelse(dfr$strand == "+", dfr$breakpoint, dfr$position))
  gr <- GRanges(seqnames = dfr$chr,
                ranges = ranges,
                strand = dfr$strand)
  
  GTSP.list <- strsplit(dfr$sampleName, split="-")
  mcols <- data.frame("count" = dfr$count,
                      "siteID" = dfr$siteID,
                      "sampleName" = dfr$sampleName,
                      "sampleID" = dfr$sampleID,
                      "GTSP" = sapply(1:length(gr), function(i){GTSP.list[[i]][1]}),
                      "refGenome" = dfr$refGenome,
                      "sex" = dfr$gender,
                      "miseqid" = dfr$miseqid)
  mcols(gr) <- mcols
  gr
}