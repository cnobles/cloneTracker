standardize_intsites <- function(sites.unstandardized, maxgap=5L, 
                               grouping=NULL, keep.mcols=FALSE, ...){
  sites.unstandardized <- sort(sites.unstandardized)
  if(is.null(grouping)){
    sites.gp <- list(sites.unstandardized)
  }else if(grouping %in% names(mcols(sites.unstandardized))){
    groups <- mcols(sites.unstandardized)[
      grep(grouping, names(mcols(sites.unstandardized)))]
    sites.unstandardized$groups <- groups[,1]
    sites.gp <- split(sites.unstandardized, sites.unstandardized$groups)
  }else{
    stop("Grouping partitioning failed. Make sure grouping is either NULL or 
         refering to the correct column in GRanges object.")
  }
  
  sites.std.grl <- GRangesList(lapply(1:length(sites.gp), function(i){
    sites <- sites.gp[[i]]
    sites$calledStart <- ifelse(strand(sites) == "+", start(sites), end(sites))
    sites$breakpoint <- ifelse(strand(sites) == "+", end(sites), start(sites))
    
    #Manipulate ranges to only have unique starts (intSites), remove breakpoint and abundance info.
    sites <- flank(sites, width = -1, start = TRUE)
    sites.rd <- reduce(sites, min.gapwidth = 0L, with.revmap = TRUE)
    sites.rd$calledStart <- ifelse(strand(sites.rd) == "+", start(sites.rd), end(sites.rd))
    sites.rd$freq <- sapply(sites.rd$revmap, length)
    #overlaps <- findOverlaps(sites.rd, maxgap = maxgap, select = "all", 
    #                         ignoreSelf = FALSE, ignoreRedundant = FALSE)
    #edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
    #clusters <- clusters(graph.edgelist(edgelist, directed = FALSE))
    clusters <- clusters(graphOverlaps(sites.rd, maxgap = maxgap))
    
    if(length(clusters$membership) > 0){
      sites.rd$clusterID <- paste0(i, ":", clusters$membership)
      sites.rd <- split(sites.rd, sites.rd$clusterID)
      sites.grl <- GRangesList(lapply(1:length(sites.rd), function(j){
        clus <- sites.rd[[j]]
        
        #At this point, position and frequency could be used to look at distribution
        #Currently, the logic uses the same as hiReadsProcessor::clusterSites(), which
        #favors the highest frequency (likely most common). I've added a mean calculation
        #for sites that are clustered together spanning greater than the windowSize, this
        #should be looked at more closely with test data sets. Lastly, in the case of equivalent
        #frequency and within the windowSize, the upstream position is prefered, as with previous
        #functions. This clusterStart calling logic should be looked at more closely, with test
        #data sets.
        
        top.freq <- clus[clus$freq == max(clus$freq),]
        if(length(top.freq) == 1){
          clus.position <- unique(top.freq$calledStart)
        }else if((range(top.freq$calledStart)[2] - range(top.freq$calledStart)[1]) > maxgap){
          message("Possible bimodal distribution of intSites in cluster ", i, ":", j, ".")
          clus.position <- as.integer(median(top.freq$calledStart))
        }else{
          clus.position <- min(top.freq$calledStart)
        }
        clus$clusterStart <- clus.position
        clus
      }))
    }else{
      message("No sites within maxgap distance, no clustering needed")
      sites.rd$clusterID <- paste0(i, ":", seq(1:length(sites.rd)))
      sites.rd$clusterStart <- sites.rd$calledStart
      sites.grl <- GRangesList(sites.rd)
    }
    
    sites.clus <- unlist(sites.grl)
    
    sites <- sites[unlist(sites.clus$revmap)]
    sites$clusterStart <- as.integer(Rle(values = sites.clus$clusterStart,
                                            lengths = sites.clus$freq))
    
    ranges <- IRanges(start = ifelse(strand(sites) == "+", 
                                     sites$clusterStart, sites$breakpoint),
                      end = ifelse(strand(sites) == "+",
                                   sites$breakpoint, sites$clusterStart))
    sites.std <- GRanges(seqnames = seqnames(sites),
                         ranges = ranges,
                         strand = strand(sites),
                         seqinfo = seqinfo(sites))
    
    if(keep.mcols){
      mcols(sites.std) <- mcols(sites[unlist(sites.clus$revmap)])
      sites.std$calledStart <- as.integer(sites$calledStart)
      sites.std$breakpoint <- as.integer(sites$breakpoint)
      sites.std$clusterStart <- as.integer(Rle(values = sites.clus$clusterStart,
                                               lengths = sites.clus$freq))
      sites.std$clusterID <- as.character(Rle(values = sites.clus$clusterID,
                                              lengths = sites.clus$freq))
    }
    
    sites.std
  }))
  sites.standardized <- unlist(sites.std.grl)
  if(!is.null(grouping)){sites.standardized$groups <- NULL}
  sites.standardized
}
