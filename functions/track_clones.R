#Return all clones present in more than 1 set of data from a list of sites
#Function makes every pairwise comparison posible from the list given
#From reading clusterSites function, which is findOverlaps based,
#I found that you don't have to specify a subject if query = subject
#and that adding the variable ignore.self (or something like that)
#which could reduce the length and demand of cloneTracker towards the end

track_clones <- function(sites, maxgap=5L, track.origin=TRUE, ...){
  grl.sites <- sites
  
  if(class(sites) == "list"){grl.sites <- GRangesList(sites)}
  
  if(track.origin == TRUE){
    list.sites <- lapply(1:length(sites), function(i){sites[[i]]})
    list.sites <- lapply(1:length(list.sites), function(i){
      list.sites[[i]]$origin <- names(sites[i])
      return(list.sites[[i]])
    })
    grl.sites <- GRangesList(list.sites)
  }
  
  condensed.sites <- unlist(grl.sites, use.names = FALSE)
  names(condensed.sites) <- 1:length(condensed.sites)
  
  #rather try findOverlaps(condensed.sites, ignoreSelf=TRUE, ignoreRedundant = FALSE, select = "all", maxgap = maxgap)    
  overlaps <- findOverlaps(condensed.sites, condensed.sites, maxgap = maxgap)
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
