#Normalize the position of integration sites between repliates and samples
normalize_intsite_positions <- function(sites, gap = 5L){
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites <- flank(sites, -1, start = TRUE)
  
  graph.gap <- graphOverlaps(sites, gap = gap)
  clusters.gap <- clusters(graph.gap)
  clusters_to_normalize <- grep(TRUE, clusters.gap$csize > 1)
  
  sites$pos.clus <- clusters.gap$membership

  norm.clus.dfl <- lapply(clusters_to_normalize, function(pos.clus){
    obs <- grep(TRUE, sites$pos.clus == pos.clus)
    dfr <- data.frame("strand" = strand(sites[obs]),
                      "called.pos" = sites[obs]$called.pos,
                      "pos.clus" = sites[obs]$pos.clus,
                      stringsAsFactors = FALSE)
    dfr
  })
  
  #Generate a key of the clusters and normalized positions
  norm.key <- as.data.frame(bind_rows(lapply(norm.clus.dfl, function(pos.dfr){
    pos.freq <- as.data.frame(table(pos.dfr$called.pos))
    top.freq <- as.numeric(as.character(
      pos.freq[pos.freq$Freq == max(pos.freq$Freq), "Var1"]
    ))
    norm.pos <- ifelse(unique(pos.dfr$strand) == "+", min(top.freq), max(top.freq))
    key.row <- data.frame("norm.pos" = norm.pos,
                          "pos.clus" = unique(pos.dfr$pos.clus),
                          stringsAsFactors = FALSE)
    key.row
  })))
  unnorm.sites <- grep(TRUE, !sites$pos.clus %in% clusters_to_normalize)
  unnorm.key <- data.frame("norm.pos" = sites[unnorm.sites]$called.pos,
                           "pos.clus" = sites[unnorm.sites]$pos.clus,
                           stringsAsFactors = FALSE)
  clus.pos.key <- as.data.frame(bind_rows(norm.key, unnorm.key))
  row.names(clus.pos.key) <- clus.pos.key$pos.clus
  sites$norm.pos <- clus.pos.key[as.character(sites$pos.clus), "norm.pos"]
  
  ranges <- IRanges(start = ifelse(strand(sites) == "+", 
                                   sites$norm.pos, 
                                   sites$called.bp),
                    end = ifelse(strand(sites) == "+", 
                                 sites$called.bp, 
                                 sites$norm.pos))
  normalized.sites <- GRanges(seqnames = seqnames(sites),
                              ranges = ranges,
                              strand = strand(sites),
                              seqinfo = seqinfo(sites))
  mcols(normalized.sites) <- mcols(sites)
  normalized.sites
}