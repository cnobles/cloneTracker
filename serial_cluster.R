#Determine cluster membership for multiple clustering gaps
serial_cluster <- function(sites, maxgaps = c(0L, 1L, 2L)){
  sites$order <- seq(1:length(sites))
  mcols <- mcols(sites)
  sites.fl <- flank(sites, -1, start = TRUE)
  serial.dfr <- bind_cols(lapply(maxgaps, function(maxgap){
    if(maxgap == 0L){
      graph <- graphOverlaps(sites.fl, maxgap = maxgap)
      membership <- clusters(graph)$membership
    }else{
      sites.rd <- reduce(sites.fl, min.gapwidth = 0L, with.revmap = TRUE)
      graph <- graphOverlaps(sites.rd, maxgap = maxgap)
      clusters <- clusters(graph)
      counts <- sapply(sites.rd$revmap, length)
      sites.reorder <- sites[unlist(sites.rd$revmap)]
      sites.reorder$membership <- as.numeric(Rle(clusters$membership, counts))
      reorder.dfr <- as.data.frame(mcols(sites.reorder))
      reorder.dfr <- arrange(reorder.dfr, order)
      membership <- reorder.dfr$membership
    }
    membership <- data.frame(membership)
    membership
  }))
  names(serial.dfr) <- sapply(maxgaps, function(x){paste0("clus.", x, "L")})
  mcols <- bind_cols(as.data.frame(mcols), serial.dfr)
  mcols$order <- NULL
  mcols(sites) <- mcols
  sites
}