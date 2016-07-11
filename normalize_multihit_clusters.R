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
      "multihitID" = as.character(gr$multihitID),
      "clusID" = gr$pos.clus,
      stringsAsFactors = FALSE
    )) 
    
    split.key <- split(key, key$clusID)
    edgelist <- matrix(c(
      Rle(
        values = sapply(split.key, function(x) x$multihitID[1]),
        lengths = sapply(split.key, nrow)
      ),
      unlist(sapply(split.key, function(x) x$multihitID))),
      ncol = 2
    )
    graph <- graph.edgelist(edgelist, directed = FALSE)
    std.multi.id <- clusters(graph)$membership
    key$stdmultihitID <- std.multi.id[key$multihitID]
    key$clusID <- NULL
    key <- unique(key)
    
    std.key <- data.frame(
      row.names = key$multihitID,
      "stdmultihitID" = key$stdmultihitID
    )
    
    gr$stdmultihitID <- std.key[as.character(gr$multihitID), "stdmultihitID"]    
    gr
  })
  
  std.multi.gr <- do.call(c, lapply(1:length(std.multi.gr), function(i){std.multi.gr[[i]]}))
  if(!is.null(grouping)){std.multi.gr$groups <- NULL}
  std.multi.gr
}
