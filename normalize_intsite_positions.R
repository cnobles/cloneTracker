#Normalize the position of integration sites between repliates and samples
normalize_intsite_positions <- function(sites, gap = 5L){
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites <- flank(sites, -1, start = TRUE)
  
  graph.gap <- graphOverlaps(sites, gap = gap)
  clusters.gap <- clusters(graph.gap)
  key.gap <- data.frame("clus" = seq(1:clusters.gap$no),
                        "size" = clusters.gap$csize)
  clusters_to_normalize <- grep(TRUE, clusters.gap$csize > 1)
  
  sites$clus <- clusters.gap$membership
  norm.sites <- sites[grep(TRUE, sites$clus %in% clusters_to_normalize)]
  
  norm.sites <- unlist(GRangesList(lapply(unique(norm.sites$clus), function(clus){
    gr <- norm.sites[norm.sites$clus == clus]
    graph.0L <- graphOverlaps(gr, gap = 0L)
    clusters.0L <- clusters(graph.0L)
    gr$freq <- as.numeric(Rle(clusters.0L$csize, clusters.0L$csize))
    gr$cliq <- clusters.0L$membership
    max.cliq <- largest.cliques(graph.0L)
    positions <- start(gr[unlist(max.cliq)])
    gr$clus.pos <- round(median(positions))
    gr
  })))
  
  ranges <- IRanges(start = ifelse(strand(norm.sites) == "+", 
                                   norm.sites$clus.pos, 
                                   norm.sites$called.bp),
                    end = ifelse(strand(norm.sites) == "+", 
                                 norm.sites$called.bp, 
                                 norm.sites$clus.pos))
  normalized.sites <- GRanges(seqnames = seqnames(norm.sites),
                       ranges = ranges,
                       strand = strand(norm.sites),
                       seqinfo = seqinfo(norm.sites))
  normalized.sites$called.pos <- norm.sites$called.pos
  normalized.sites$called.bp <- norm.sites$called.bp
  normalized.sites$clus <- norm.sites$clus
  normalized.sites <- c(normalized.sites, 
                        sites[!sites$clus %in% clusters_to_normalize])
  normalized.sites <- sort(normalized.sites)
  normalized.sites
}