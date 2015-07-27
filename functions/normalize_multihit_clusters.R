normalize_multihit_clutsers <- function(multihit.gr){
  #Multihits must be standardized and have clusterID info
  key <- data.frame("multihitID" = multi.gr$multihitID,
                    "clusterID" = as.character(multi.gr$clusterID)) 
  
  clusterID_to_position <- data.frame("clusterID" = unique(key$clusterID), 
                                      "position" = seq(1:length(unique(key$clusterID))))
  
  std.clusID <- merge(key, clusterID_to_position, by.x="clusterID", by.y="clusterID")
  
  multi.art.gr <- GRanges(seqnames = Rle("art_chr", nrow(key)),
                          ranges = IRanges(start = key$position, width = 1),
                          strand = "*",
                          multihitID = key$multihitID,
                          clusterID = key$clusterID)
  
  multi.art.gl <- split(multi.art.gr, multi.art.gr$multihitID)
  overlaps <- findOverlaps(multi.art.gl, minoverlap = 0L, select = "all", 
                           ignoreSelf = FALSE, ignoreRedundant = FALSE)
  edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
  clusters <- clusters(graph.edgelist(edgelist, directed = FALSE))
  
  clus.dfr <- data.frame("multihitID" = names(multi.art.gl),
                         "multihitID2" = clusters$membership)
  
  key <- merge(key, clus.dfr, by.x = "multihitID", by.y = "multihitID")
  multi.list <- split(multi.gr, multi.gr$multihitID)
  multi.list <- lapply(multi.list, function(x){
    multihitID <- unique(x$multihitID)
    multihitID2 <- key[key$multihitID == multihitID, 4]
    x$multihitID2 <- multihitID2
    x
  })
  multi.gr <- do.call(c, lapply(1:length(multi.list), function(i){multi.list[[i]]}))
  multi.gr
}