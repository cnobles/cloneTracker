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


standardize_by_clustering <- function(unstandardized.sites, standardize_breakpoints = TRUE, 
                                      get.analysis.data = TRUE, calc.cluster.stats = TRUE){
  #Build data.frame with called info and cluster information
  raw.sites <- unstandardized.sites
  raw.positions <- flank(granges(raw.sites), -1, start = TRUE)
  raw.breakpoints <- flank(granges(raw.sites), -1, start = FALSE)
  raw.positions <- serial_cluster(raw.positions, maxgaps = c(1L, 5L))
  raw.breakpoints <- serial_cluster(raw.breakpoints, maxgaps = c(1L))
  clus.dfr <- data.frame(
    "seqnames" = seqnames(raw.sites),
    "strand" = strand(raw.sites),
    "called.pos" = start(raw.positions),
    "called.bp" = end(raw.breakpoints),
    "pos.clus.5L" = raw.positions$clus.5L,
    "pos.clus.1L" = raw.positions$clus.1L,
    "bp.clus.1L" = raw.breakpoints$clus.1L
    )
  
  #Using the cluster membership for both ends of the reads, determine the
  #standardized membership for both ends
  #I should compare results, but I may not need to split this data up
  bp.pos.clus <- lapply(unique(clus.dfr$bp.clus.1L), function(bp){
    unique(clus.dfr[clus.dfr$bp.clus.1L == bp,]$pos.clus.1L)
  })
  pos.bp.clus <- lapply(unique(clus.dfr$pos.clus.1L), function(pos){
    unique(clus.dfr[clus.dfr$pos.clus.1L == pos,]$bp.clus.1L)
  })
  pos.el <- as.matrix(bind_rows(lapply(unique(clus.dfr$pos.clus.1L), function(pos){
    bps <- pos.bp.clus[[pos]]
    subject <- unique(unlist(bp.pos.clus[bps]))
    dfr <- data.frame("query.pos" = pos, "subject.pos" = subject)
    dfr
  })), ncol = 2)
  pos.graph <- graph.edgelist(pos.el)
  adj.pos.clus <- clusters(pos.graph)$membership
  clus.dfr$adj.pos.clus <- adj.pos.clus[clus.dfr$pos.clus.1L]

  #Determine the standardized position and breakpoints independently
  #Positions
  pos.dfr <- clus.dfr[,c("seqnames", "strand","called.pos", "adj.pos.clus", "bp.clus.1L")]
  pos.list <- split(pos.dfr, pos.dfr$adj.pos.clus)
  std.positions <- as.data.frame(bind_rows(lapply(pos.list, function(dfr){
    dfr <- distinct(dfr)
    pos.freq <- as.data.frame(table(dfr$called.pos))
    top.freq <- as.numeric(as.character(pos.freq[pos.freq$Freq == max(pos.freq$Freq), "Var1"]))
    std.pos <- ifelse(unique(dfr$strand) == "+", min(top.freq), max(top.freq))
    adj.pos.clus.dfr <- data.frame("adj.pos.clus" = unique(dfr$adj.pos.clus), 
                                   "std.pos" = as.integer(std.pos))
    adj.pos.clus.dfr
  })))
  row.names(std.positions) <- std.positions$adj.pos.clus
  clus.dfr$std.pos <- std.positions[clus.dfr$adj.pos.clus, "std.pos"]
  
  #Breakpoints (optional currently, only recommended on raw ranges or within a replicate, not across replicates or samples)
  if(standardize_breakpoints){
    bp.dfr <- clus.dfr[,c("seqnames", "strand","called.bp", "adj.pos.clus", "bp.clus.1L")]
    bp.list <- split(bp.dfr, bp.dfr$bp.clus.1L)
    std.breakpoints <- as.data.frame(bind_rows(lapply(bp.list, function(dfr){
      dfr <- distinct(dfr)
      bp.freq <- as.data.frame(table(dfr$called.bp))
      top.freq <- as.numeric(as.character(bp.freq[bp.freq$Freq == max(bp.freq$Freq), "Var1"]))
      std.bp <- ifelse(unique(dfr$strand) == "+", max(top.freq), min(top.freq))
      adj.bp.clus.dfr <- data.frame("adj.bp.clus" = unique(dfr$bp.clus.1L), 
                                    "std.bp" = as.integer(std.bp))
      adj.bp.clus.dfr
    })))
    row.names(std.breakpoints) <- std.breakpoints$adj.bp.clus
    clus.dfr$std.bp <- std.breakpoints[clus.dfr$bp.clus.1L, "std.bp"]
  }else{
    clus.dfr$std.bp <- clus.dfr$called.bp
  }
  
  #Rebuild Granges of standardized sites
  std.ranges <- IRanges(start = ifelse(clus.dfr$strand == "+", clus.dfr$std.pos, clus.dfr$std.bp),
                        end = ifelse(clus.dfr$strand == "+", clus.dfr$std.bp, clus.dfr$std.pos))
  std.sites <- GRanges(seqnames = clus.dfr$seqnames,
                       ranges = std.ranges,
                       strand = clus.dfr$strand,
                       seqinfo = seqinfo(raw.sites))
  mcols(std.sites) <- mcols(raw.sites)
  
  if(get.analysis.data){
    std.sites$called.pos <- clus.dfr$called.pos
    std.sites$called.bp <- clus.dfr$called.bp
    std.sites$pos.clus.ori <- clus.dfr$pos.clus.1L
    std.sites$pos.clus.adj <- clus.dfr$adj.pos.clus
    std.sites$bp.clus <- clus.dfr$bp.clus.1L
  }
  
  if(calc.cluster.stats){
    clus.list <- split(clus.dfr, clus.dfr$adj.pos.clus)
    clus.stats <- as.data.frame(bind_rows(lapply(clus.list, function(clus){
      bp.range.diff <- diff(range(clus$called.bp))
      bp.ori.count <- length(unique(clus$called.bp))
      bp.clus.count <- length(unique(clus$bp.clus.1L))
      bp.mean.diff <- ifelse(nrow(clus) > 1, mean(diff(clus$called.bp)), 0)
      stats.dfr <- data.frame(bp.range.diff, bp.ori.count, 
                              bp.clus.count, bp.mean.diff)
      stats.dfr
    })))
    row.names(clus.stats) <- as.numeric(names(clus.list))
    
    #append the data to std.sites
    std.sites$bp.range.diff <- clus.stats[clus.dfr$adj.pos.clus,"bp.range.diff"]
    std.sites$bp.ori.count <- clus.stats[clus.dfr$adj.pos.clus, "bp.ori.count"]
    std.sites$bp.clus.count <- clus.stats[clus.dfr$adj.pos.clus, "bp.clus.count"]
    std.sites$bp.mean.diff <- clus.stats[clus.dfr$adj.pos.clus, "bp.mean.diff"]
  }
  
  return(std.sites)
}





















































