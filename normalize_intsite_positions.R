#Normalize the position of integration sites between repliates and samples
normalize_intsite_positions <- function(sites, gap = 5L){
  sites$order <- 1:length(sites)
  mcols <- mcols(sites)
  mcols(sites) <- NULL
  sites$order <- 1:length(sites)
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  mcols$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  mcols$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites <- flank(sites, -1, start = TRUE)
  
  graph.gap <- graphOverlaps(sites, gap = gap)
  clusters.gap <- clusters(graph.gap)
  clusters_to_normalize <- grep(TRUE, clusters.gap$csize > 1)
  
  sites$clus <- clusters.gap$membership
  mcols$clus <- clusters.gap$membership
  
  unnorm.sites <- sites[grep(TRUE, !sites$clus %in% clusters_to_normalize)]
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
  
  normalized.sites$order <- norm.sites$order
  unnorm.sites$called.pos <- NULL
  unnorm.sites$called.bp <- NULL
  unnorm.sites$clus <- NULL
  normalized.sites <- c(normalized.sites, unnorm.sites)
  
  mcols <- as.data.frame(mcols, row.names = mcols$order)[normalized.sites$order,]
  mcols(normalized.sites) <- mcols
  normalized.sites <- unlist(GRangesList(sapply(1:length(normalized.sites),
                                                function(i){
                                                  normalized.sites[normalized.sites$order == i]
                                                })))
  normalized.sites$order <- NULL
  normalized.sites
}