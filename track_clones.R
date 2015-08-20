#Return all clones present in more than 1 set of data from a list of sites

track_clones <- function(sites.list, gap=5L, track.origin=TRUE, ...){
  grl.sites <- sites.list
  
  if(class(sites.list) == "list"){grl.sites <- GRangesList(sites.list)}
  
  if(track.origin == TRUE){
    list.sites <- lapply(1:length(sites.list), function(i){sites.list[[i]]})
    list.sites <- lapply(1:length(list.sites), function(i){
      list.sites[[i]]$origin <- names(sites.list[i])
      return(list.sites[[i]])
    })
    grl.sites <- GRangesList(list.sites)
  }
  
  condensed.sites <- unlist(grl.sites, use.names = FALSE)
  names(condensed.sites) <- 1:length(condensed.sites)
  
  overlaps <- findOverlaps(condensed.sites, maxgap = gap, ignoreSelf = FALSE, ignoreRedundant = FALSE)
  edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
  
  clusters <- clusters(graph.edgelist(edgelist, directed = FALSE))
  clusters.names <- split(names(condensed.sites), clusters$membership)
  clusters.lengths <- data.frame(
    id = c(1:length(clusters.names)),
    length = sapply(clusters.names, function(x){
      length(x)
    }))
  
  clusters.true <- clusters.names[
    clusters.lengths[clusters.lengths$length > 1,"id"]]
  
  clustered.sites <- GRangesList(lapply(clusters.true, function(x){
    unname(condensed.sites[x])
  }))
  
  return(clustered.sites)
}
