normalize_multihit_clutsers <- function(multihits.gr, gap = 5L, grouping = NULL){
  #Multihits must be standardized and have clusterID info
  if(is.null(grouping)){
    multihits.gp <- GRangesList(multihits.gr)
  }else if(grouping %in% names(mcols(multihits.gr))){
    groups <- mcols(multihits.gr)[
      grep(grouping, names(mcols(multihits.gr)))]
    multihits.gr$groups <- groups[,1]
    multihits.gp <- split(multihits.gr, multihits.gr$groups)
  }else{
    stop("Grouping partitioning failed. Make sure grouping is either NULL or 
         refering to the correct column in GRanges object.")
  }
  
  std.multi.gr <- lapply(multihits.gp, function(gr){
    key <- unique(data.frame(
      "multihitID" = gr$multihitID,
      "clusID" = gr$pos.clus
    )) 
  
    #clusterID_to_position <- data.frame("clusterID" = unique(key$clusterID), 
    #                                    "position" = seq(1:length(unique(key$clusterID))))
  
    #key <- merge(key, clusterID_to_position, by.x="clusterID", by.y="clusterID")
    #key <- unique(key)
  
    split.key <- split(key, key$clusID)
    edgelist <- matrix(c(
      as.numeric(Rle(
        values = sapply(split.key, function(x) x$multihitID[1]),
        lengths = sapply(split.key, nrow)
      )),
      unlist(sapply(split.key, function(x) x$multihitID))),
      ncol = 2
    )
    graph <- graph.edgelist(edgelist, directed = FALSE)
    key$stdmultiID <- clusters(graph)[[1]]
    
    
    multi.art.gr <- GRanges(seqnames = Rle("art_chr", nrow(key)),
                            ranges = IRanges(start = key$position, width = 1),
                            strand = Rle("*", nrow(key)),
                            multihitID = key$multihitID,
                            clusterID = key$clusterID)
    graph <- graphOverlaps(multi.art.gr, gap = 0L)
    key$stdclusID <- clusters(graph)$membership
    
    el <- matrix( c(multi.art.gr$multihitID, stdmultihitID), ncol = 2 )
    graph.el <- graph.edgelist(el, directed = FALSE)
    clus <- clusters(graph.el)
    
    
    #multi.art.gl <- split(multi.art.gr, multi.art.gr$multihitID)
    #overlaps <- findOverlaps(multi.art.gl, minoverlap = 0L, select = "all", 
    #                       ignoreSelf = FALSE, ignoreRedundant = FALSE)
    #edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
    #clusters <- clusters(graph.edgelist(edgelist, directed = FALSE))
  
    #clus.dfr <- data.frame("multihitID" = names(multi.art.gl),
    #                     "stdmultihitID" = clusters$membership)
    clus.dfr <- data.frame("multihitID" = multi.art.gr$multihitID,
                           "stdmultihitID" = stdmultihitID)
    
    key <- merge(key, clus.dfr, by.x = "multihitID", by.y = "multihitID")
    multi.list <- split(gr, gr$multihitID)
    multi.list <- lapply(multi.list, function(x){
      multihitID <- unique(x$multihitID)
      stdmultihitID <- key[key$multihitID == multihitID, 4]
      x$stdmultihitID <- stdmultihitID
      x
    })
    gr <- do.call(c, lapply(1:length(multi.list), function(i){multi.list[[i]]}))
    gr
  })
  
  std.multi.gr <- do.call(c, lapply(1:length(std.multi.gr), function(i){std.multi.gr[[i]]}))
  if(!is.null(grouping)){std.multi.gr$groups <- NULL}
  std.multi.gr
}
