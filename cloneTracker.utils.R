#Using a keep_cols list, remove unwanted metadata from GRanges
condense_metadata <- function(sites, keep_cols){
  tot_metadata <- names(mcols(sites))
  keep_positions <- match(keep_cols, tot_metadata)
  cleaned_sites <- sites[, keep_positions]
  return(cleaned_sites)
}


#Clean a GRange object to only have the range information
clean_sites <- function(sites){
  cleaned_sites <- granges(sites)
  names(cleaned_sites) <- NULL
  cleaned_sites
}


#Generate position ID (posid) given intsite parameters for class::GRange
generate_posID <- function(sites=NULL, seqnames=NULL, strand=NULL, start=NULL, end=NULL, ...){
  if(length(sites) != 0){
    if(class(sites) == "GRanges"){
      chr <- as.character(seqnames(sites))
      strand <- as.vector(strand(sites))
      pos <- ifelse(strand == "+", start(sites), end(sites))
      posID <- paste0(chr, strand, pos)
    }else{
      message("Sites provided not a GRanges object, please use alternative inputs.")
      stop()
    }
  }else{
    if(length(seqnames) != 0 & length(strand) != 0 & length(start) != 0 & length(end) != 0){
      chr <- as.character(seqnames)
      strand <- as.vector(strand)
      start <- as.integer(start)
      end <- as.integer(end)
      sites.df <- data.frame(chr, strand, start, end)
      sites.df$pos <- ifelse(strand == "+", sites.df$start, sites.df$end)
      posID <- paste0(sites.df$chr, sites.df$strand, sites.df$pos)
    }else{
      message("Please supply seqnames, strand, start, and end info.")
      stop()
    }}
  return(posID)
}


#Generate graph of overlaping (determined by maxgap size) GRange observations
graphOverlaps <- function(sites, gap){
  overlaps <- findOverlaps(sites, maxgap = gap)
  edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
  graph <- graph.edgelist(edgelist, directed = FALSE)
  graph
}


#Find all sites/clones using a list of position ID's (posid)
find_sites <- function(sites, posID){
  request <- do.call(c, lapply(1:length(posID), function(i){
    sites[sites$posID == posID[[i]],]
  }))
}



#Developmental functions, not currently used
#Difference between this function and flank(x, -1, start=TRUE) is neglegible
.intSiteCollapse <- function(sites, use.names=FALSE){
  if(use.names == TRUE){if(length(names(sites)) == 0){
    message("Sites do not have names. Changing use.names = FALSE")
    use.names <- FALSE
  }}
  
  strand <- as.vector(strand(sites))
  adj.start <- ifelse(strand == "+", start(sites), end(sites))
  adj.width <- rep(1, length(sites))
  adj.ranges <- IRanges(start = adj.start, width = adj.width)
  adj.sites <- GRanges(seqnames = seqnames(sites), 
                       ranges = adj.ranges, strand = strand(sites))
  
  mcols(adj.sites) <- mcols(sites)
  if(use.names == TRUE){names(adj.sites) <- names(sites)}
  
  adj.sites
}


.intSiteCluster <- function(sites, windowSize=5L, grouping=1){
  
  clusters <- clusterSites(
    posID = paste0(seqnames(sites), "_", strand(sites)),
    value = start(sites),
    grouping = grouping,
    windowSize = windowSize)
  
  seqnames_strand <- strsplit(clusters$posID, split="_")
  
  sites.df <- data.frame(
    seqnames = sapply(1:length(seqnames_strand), function(i){seqnames_strand[[i]][1]}),
    strand = sapply(1:length(seqnames_strand), function(i){seqnames_strand[[i]][2]}),
    start = clusters[,1])
  
  sites.df <- distinct(sites.df)
  
  sites.clustered <- GRanges(
    seqnames = as.character(sites.df$seqnames),
    ranges = IRanges(
      start = sites.df$start, 
      width = rep(1, nrow(sites.df))),
    strand = as.character(sites.df$strand))
  
  return(sites.clustered)
}

.standardize_by_overlaps <- function(sites.unstandardized, gap=5L, 
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
    clusters <- clusters(graphOverlaps(sites.rd, gap = gap))
    
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
        }else if((range(top.freq$calledStart)[2] - range(top.freq$calledStart)[1]) > gap){
          message("Possible bimodal distribution of intSites in cluster ", i, ":", j, ".")
          clus.position <- as.integer(median(top.freq$calledStart))
        }else{
          clus.position <- min(top.freq$calledStart)
        }
        clus$clusterStart <- clus.position
        clus
      }))
    }else{
      message("No sites within gap distance, no clustering needed")
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


#Additional functionality over .standardizeSites.
.standardize_by_clusSites <- function(sites, windowSize=5L, grouping=1, condense=FALSE, 
                                pcrBreakpoints=FALSE, returnAbund=FALSE, 
                                abundMethod="fragLen", keepMcols=FALSE, ...){
  
  if(keepMcols == TRUE & length(mcols(sites)) == 0){keepMcols <- FALSE}
  
  if(returnAbund){
    if(abundMethod == "fragLen"){
      Abund <- FALSE
    }else if(abundMethod == "estAbund"){
      Abund <- TRUE
    }else{
      stop("Must choose either fragLen or estAbund for abundMethod.")
    }
  }else{
    Abund <- FALSE
  }
  
  sites <- sort(sites)
  mcols <- mcols(sites)
  mcols(sites) <- NULL
  
  #Add required columns to pass into clusterSites as a psl.rd
  sites$Position <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$Break <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites$score <- rep(95, length(sites))
  sites$qEnd <- width(sites)
  
  #Positions clustered by 5L window and best position is chosen for cluster
  sites.standardized <- clusterSites(
    psl.rd = sites,
    weight = rep(1, length(sites)),
    sonicAbund = Abund)
  
  start(sites.standardized) <- ifelse(strand(sites.standardized) == "+", 
                                      sites.standardized$clusteredPosition, 
                                      sites.standardized$Break)
  end(sites.standardized) <- ifelse(strand(sites.standardized) == "-",
                                    sites.standardized$clusteredPosition, 
                                    sites.standardized$Break)
  
  sites.standardized$Position <- NULL
  sites.standardized$Break <- NULL
  sites.standardized$score <- NULL
  sites.standardized$qEnd <- NULL
  sites.standardized$clusteredPosition <- NULL
  if(condense){
    sites.standardized$clonecount <- NULL
    sites.standardized$clusterTopHit <- NULL
  }
  sites.standardized <- sort(sites.standardized)
  
  if(returnAbund){
    if(abundMethod == "fragLen"){
      estAbund.uniqueFragLen <- function(location, fragLen, replicate=NULL){
        if(is.null(replicate)){replicate <- 1}  #Need for downstream workflow
        dfr <- data.frame(location = location, fragLen = fragLen, 
                          replicate = replicate)
        dfr_dist <- distinct(dfr)
        site_list <- split(dfr_dist, dfr_dist$location)
        theta <- sapply(site_list, function(x){nrow(x)})
        theta <- theta[unique(dfr$location)]
        list(theta=theta)
      }
      
      reps <- which("replicate" == names(mcols))
      if(length(reps) != 0){
        replicate <- mcols$replicate
      }else{
        replicate <- rep(1, length(sites.standardized))
      }
      
      posid <- generate_posID(sites.standardized)
      dfr <- data.frame("ID" = posid,
                        "fragLength" = width(sites.standardized),
                        "replicate" = replicate)
      
      if(length(unique(dfr$replicate)) == 1){
        estimatedAbundances <- estAbund.uniqueFragLen(dfr$ID, dfr$fragLength)
      }else{
        estimatedAbundances <- estAbund.uniqueFragLen(dfr$ID, dfr$fragLength, dfr$replicate)
      }
      
      sites.standardized$estAbund <- round(estimatedAbundances$theta)
      sites.standardized$estAbundProp <- sites.standardized$estAbund/sum(sites.standardized$estAbund)
      sites.standardized$estAbundRank <- rank(-1*sites.standardized$estAbundProp, ties.method="max")
      
      rm(estimatedAbundances, dfr, posid, replicate)
    }
  }
  
  if(keepMcols == FALSE){mcols <- NULL}
  if(returnAbund){
    if(is.null(mcols)){
      mcols <- mcols(sites.standardized)
    }else{
      mcols <- cbind(mcols, mcols(sites.standardized))
    }
  }
  mcols(sites.standardized) <- mcols
  
  if(condense){
    sites.reduced <- intSiteCollapse(sites.standardized)
    sites.reduced <- unlist(reduce(sites.reduced, with.revmap=TRUE))
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    sites.condensed <- sites.standardized[unlist(sites.reduced$revmap)]
    sites.condensed <- split(sites.condensed, Rle(values = seq(length(sites.reduced)), 
                                                  lengths = sites.reduced$counts))
    
    #Try to change this into an integrer list rather than character.string
    if(pcrBreakpoints == TRUE){
      comma.paste <- function(...){paste(..., sep=",")}
      breakpoints <- lapply(sites.condensed, width)
      breakpoints <- lapply(1:length(breakpoints), function(i){
        do.call(comma.paste, lapply(1:length(breakpoints[[i]]), 
                                    function(j){breakpoints[[i]][j]}))})
      breakpoints <- unlist(breakpoints)
    }
    
    mcols.condensed <- lapply(1:length(sites.condensed), function(i){
      mcols.i <- mcols(sites.condensed[[i]][1])})
    mcols.condensed <- do.call(rbind, lapply(1:length(mcols.condensed), function(i){
      mcols.condensed[[i]]}))
    
    sites.condensed <- unlist(reduce(sites.condensed))
    mcols(sites.condensed) <- mcols.condensed
    
    if(pcrBreakpoints == TRUE){sites.condensed$pcrBreakpoints <- breakpoints}
  }
  
  if(condense == FALSE){
    sites.requested <- sites.standardized
  }else{
    sites.requested <- sites.condensed
  }
  
  return(sites.requested)
}

.standardizeSites <- function(unstandardizedSites){
  if( ! length(unstandardizedSites) > 0){
    return(unstandardizedSites)
  }
  #Get called start values for clustering  
  unstandardizedSites$Position <- ifelse(strand(unstandardizedSites) == "+", start(unstandardizedSites), end(unstandardizedSites))
  unstandardizedSites$Break <- ifelse(strand(unstandardizedSites) == "+", end(unstandardizedSites), start(unstandardizedSites))
  unstandardizedSites$Score <- 95
  unstandardizedSites$qEnd <- width(unstandardizedSites)
  
  #Positions clustered by 5L window and best position is chosen for cluster
  standardized <- clusterSites(
    psl.rd = unstandardizedSites,
    weight = rep(1, length(unstandardizedSites)) 
  )
  
  start(standardized) <- ifelse(strand(standardized) == "+", 
                                standardized$clusteredPosition, standardized$Break)
  end(standardized) <- ifelse(strand(standardized) == "-", 
                              standardized$clusteredPosition, standardized$Break)
  
  standardized$Position <- NULL
  standardized$Break <- NULL
  standardized$score <- NULL
  standardized$qEnd <- NULL
  standardized$clusteredPosition <- NULL
  standardized$clonecount <- NULL
  standardized$clusterTopHit <- NULL
  
  sort(standardized)
}


#A logic table simply reports on the presence or absence of a clone in a 
#data set
.generate_logic_table <- function(hits, list){
  table <- do.call(cbind, lapply(1:length(list), function(j){
    logic_col <- sapply(1:length(hits), function(i){
      logic <- findOverlaps(hits[i], list[[j]], maxgap = 5L, minoverlap = 1L,
                            type = "any")
      if(length(logic) != 0){outcome <- TRUE}else{outcome <- FALSE}
      return(outcome)
    })
    return(logic_col)
  }))
  colnames(table) <- names(list)
  rownames(table) <- hits$posid
  return(table)
}

#Generate a table which contains all the sonic abundances for each site
.generate_abund_table <- function()
  
  
#Though likely modified for specific situations, score based on logic table
#to determine best ranking or tracked clones
.score_tracking_hits <- function(i, table){
  score <- 0
  
  if(table[i,1] == TRUE){
    score <- add(1, score)
  }else{if(table[i,2] == TRUE){
    score <- add(1, score)
  }else{
    score <- add(0, score)
  }}
  
  if(table[i,3] == TRUE){
    score <- add(2, score)
  }else{if(table[i,4] == TRUE){
    score <- add(2, score)
  }else{
    score <- add(0, score)
  }}
  
  if(table[i,5] == TRUE){
    score <- add(3, score)
  }else{
    score <- add(0, score)
  }
  
  return(score)}