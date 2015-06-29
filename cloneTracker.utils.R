#Dependancies check, if required packages are not loaded, script is haulted
dependancies <- c("dplyr", "IRanges", "GenomicRanges", "igraph", "hiReadsProcessor")

sapply(dependancies, function(package){
  library(package, character.only = TRUE)})

dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(FALSE %in% dependancies_present){
  message("Load required packages. Check Unloaded_Packages for missing
          dependancies.")
  Unloaded_Packages <- data.frame(package=as.character(dependancies), 
                                  loaded=dependancies_present)
  stop()
  }else{
  remove(dependancies, dependancies_present)
  message("Required packages loaded.")
}

#Test GRange
test <- GRanges(seqnames = Rle(values = c("chr1", "chr2", "chr3"),
                               lengths = c(5, 10, 7)),
                ranges = IRanges(start = c(2,7,10,7,8,15,15,17,16,14,14,6,7,6,5,2,5,8,11,14,17,20),
                                 width = rep(30, 22)),
                strand = Rle(values = c("+","-","+"),
                             lengths = c(8,7,7)))

#Utility functions for clone tracking using integration sites as markers
#Remove width data from called intsites, returned data is used for
#clone marking

intSiteCollapse <- function(sites, use.names=FALSE){
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
 
  return(adj.sites)
}

.intSitesCluster <- function(sites, windowSize=5L, grouping=1){
  
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

#intSiteStandardize, change name at some point
.intSiteStandardize <- function(sites, windowSize=5L, grouping="1", condense=FALSE, 
                                pcr.breakpoints=FALSE, return.sonicAbund=FALSE, 
                                keep.mcols=TRUE, ...){
  sites <- sort(sites)
  true_starts <- ifelse(strand(sites) == "+", start(sites), end(sites))
  
  clusters <- clusterSites(
    posID = paste0(seqnames(sites), "_", strand(sites)),
    value = true_starts,
    grouping = grouping,
    windowSize = windowSize)
  
  seqnames_strand <- strsplit(as.character(Rle(values = clusters$posID,
                                               lengths = clusters$freq)), split="_")
  
  sites.df <- data.frame(
    seqnames = sapply(1:length(seqnames_strand), function(i){seqnames_strand[[i]][1]}),
    strand = sapply(1:length(seqnames_strand), function(i){seqnames_strand[[i]][2]}),
    start.cluster = Rle(values = clusters$clusteredValue, lengths = clusters$freq),
    start.unstand = Rle(values = clusters$value, lengths = clusters$freq),
    freq = Rle(values = clusters$freq, lengths = clusters$freq))
  
  sites.df <- arrange(sites.df, seqnames, desc(strand), start.unstand)
  sites.df$end <- ifelse(strand(sites) == "+", end(sites), start(sites))
  
  range.standardized <- IRanges(
    start = ifelse(sites.df$strand == "+", sites.df$start.cluster, sites.df$end),
    end = ifelse(sites.df$strand == "+", sites.df$end, sites.df$start.cluster))
  
  sites.standardized <- GRanges(
    seqnames = sites.df$seqnames,
    ranges = range.standardized,
    strand = sites.df$strand,
    seqinfo = seqinfo(sites))
  
  if(keep.mcols == TRUE){
    mcols <- mcols(sites)
  }else{
    mcols <- NULL
  }
  
  mcols(sites.standardized) <- mcols
  
  if(condense == TRUE){
    sites.reduced <- intSiteCollapse(sites.standardized)
    sites.reduced <- unlist(reduce(sites.reduced, with.revmap=TRUE))
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    #Not sure if I need this as the sites.standardized were sorted already, and this function would
    #just unlist to a sequence of 1:length(sites.standardized)
    #sites.condensed <- sites.standardized[unlist(sites.reduced$revmap)]
    sites.condensed <- split(sites.standardized, Rle(values = seq(length(sites.reduced)), 
                                                  lengths = sites.reduced$counts))
    
    if(pcr.breakpoints == TRUE){
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
    
    if(pcr.breakpoints == TRUE){sites.condensed$pcr.breakpoints <- breakpoints}
    
    if(return.sonicAbund == TRUE){
      posID <- generate_posid(sites.condensed)
      fragLen <- lapply(1:length(breakpoints), function(i){
        breaks <- strsplit(breakpoints[[i]], split=",")
        breaks <- as.integer(breaks[[1]])})
      sonicAbund <- do.call(rbind, lapply(1:length(sites.condensed), function(i){
        getSonicAbund(posID = posID[i], fragLen = fragLen[[i]], grouping = "1")}))
    #check getSonicAbund actually works, seems to give some weird output in estAbund for single sites
      mcols.condensed <- mcols(sites.condensed)
      mcols.condensed$posID <- generate_posid(sites.condensed)
      mcols.condensed$posID2 <- sonicAbund[,2]
      mcols.condensed$estAbund <- sonicAbund[,3]
      merge.check <- data.frame(table(mcols.condensed$posID == mcols.condensed$posID2))
      
      if(merge.check[merge.check[,1] == TRUE,"Freq"] != nrow(mcols.condensed)){
        message("Failed to merge estAbund with mcols.condensed.")
        stop()}
      
      mcols(sites.condensed) <- mcols.condensed
      
      sites.condensed$posID3 <- generate_posid(sites.condensed)
      merge.check <- data.frame(table(sites.condensed$posID == sites.condensed$posID3))
      if(merge.check[merge.check[,1] == TRUE, "Freq"] != length(sites.condensed)){
        message("Failed to correctly merge mcols with sites.condensed.")
        stop()}
      
      sites.condensed$posID <- sites.condensed$posID2 <- sites.condensed$posID3 <- NULL
    }
  }  
    
  if(condense == FALSE){
    sites.requested <- sites.standardized
  }else{
    sites.requested <- sites.condensed
  }
  
  return(sites.requested)
}




##Function to cluster integration sites into unique integration sites
#intSiteCluster <- function(sites, track.origin=TRUE, origin.column=NULL, 
#                           sonicAbund=TRUE, ...){
#  
#  if(track.origin == TRUE){
#    if(length(origin.column) == 0){
#      message("Select origin.column.")
#      stop()}}
#  
#  if(mean(width(sites)) != 1){
#    message("Collapse intSites before clustering.")
#    stop()}
#  
#  posID <- paste0(seqnames(sites), strand(sites))
#  value <- start(sites)
#  fragLen <- width(sites)
#  
#  clusters <- clusterSites(posID = posID, value = value)
  
  #Maybe just take away clusterSites and use findOverlaps and igraph (if needed)
  
  #Generate posID for clusterSites
  #from clusterSites, generate clustered posID's and pull lengths from sites 
  #using the clusterSites output as a key to link posID1 and posID2
#  }


#Return all clones present in more than 1 set of data from a list of sites
#Function makes every pairwise comparison posible from the list given
#From reading clusterSites function, which is findOverlaps based,
#I found that you don't have to specify a subject if query = subject
#and that adding the variable ignore.self (or something like that)
#which could reduce the length and demand of cloneTracker towards the end

cloneTracker <- function(sites, maxgap=5L, track.origin=TRUE, ...){
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
    
  overlaps <- findOverlaps(condensed.sites, condensed.sites, maxgap = maxgap)
  edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol = 2)
  
  clusters <- clusters(graph.edgelist(edgelist, directed = FALSE)) #Why false directed?
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


#In order to use GRange objects with hiReadsProcessor package,
#we'll need to get them in the correct format, as if they were psl files
.prepGRanges <- function(sites, origin.column=origin.column, ...){
  
  mcols <- mcols(sites)
  
  sites.psl <- data.frame(
    matches = as.numeric(width(sites)),
    misMatches = as.numeric(rep(0, length(sites))),
    repMatches = as.numeric(rep(0, length(sites))),
    nCount = as.numeric(rep(0, length(sites))),
    qNumInsert = as.numeric(rep(0, length(sites))),
    qBaseInsert = as.numeric(rep(0, length(sites))),
    tNumInsert = as.numeric(rep(0, length(sites))),
    tBaseInsert = as.numeric(rep(0, length(sites))),
    strand = as.character(strand(sites)),
    qName = as.character(mcols[, sample.name]),
    qSize = as.numeric(width(sites)),
    qStart = as.numeric(rep(0, length(sites))),
    qEnd = as.numeric(width(sites)),
    tName = as.character(seqnames(sites)),
    tSize = as.numeric(rep(NA, length(sites))),
    tStart = as.numeric(start(sites)),
    tEnd = as.numeric(end(sites)),
    blockCount = as.numeric(rep(1, length(sites))),
    blockSizes = as.character(rep(NA, length(sites))),
    qStarts = as.character(rep(NA, length(sites))),
    tStarts = as.character(rep(NA, length(sites)))
  )
  
  return(sites.psl, mcols)
}

#Find all sites/clones using a list of position ID's (posid)
find_sites <- function(sites, posid){
  request <- do.call(c, lapply(1:length(posid), function(i){
    sites[sites$posid == posid[[i]],]
  }))
}

#After clustering sites by hiReadsProcessor::clusterSites, select only
#top hits, giving the dominant site for the cluster
.collect_top_hits_for_int_sites <- function(sites){
  sites <- sites[sites$clusterTopHit == TRUE,]
  return(sites)
}

#Using a keep_cols list, remove unwanted metadata from GRanges
condense_metadata <- function(sites, keep_cols){
  tot_metadata <- names(mcols(sites))
  keep_positions <- match(keep_cols, tot_metadata)
  cleaned_sites <- sites[, keep_positions]
  return(cleaned_sites)
}

#Generate position ID (posid) given intsite parameters for class::GRange
generate_posid <- function(sites){
  chr <- as.character(seqnames(sites))
  strand <- as.vector(strand(sites))
  pos <- ifelse(strand == "+", start(sites), end(sites))
  posID <- paste0(chr, strand, pos)
  return(posID)
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