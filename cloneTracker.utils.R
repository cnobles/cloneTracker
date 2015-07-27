#Dependancies check, if required packages are not loaded, script is haulted
dependancies <- c("dplyr", "IRanges", "GenomicRanges", "igraph")

sapply(dependancies, function(package){
  library(package, character.only = TRUE)})

dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(FALSE %in% dependancies_present){
  Unloaded_Packages <- data.frame(package=as.character(dependancies), 
                                  loaded=dependancies_present)
  stop("Load required packages. Check Unloaded_Packages for missing
          dependancies.")
  }else{
  remove(dependancies, dependancies_present)
  message("Required packages loaded.")
}

#Source function scripts
function_scripts <- c("track_clones.R", 
                      "determine_abundance.R", 
                      "standardize_intsites.R", 
                      "normalize_multihit_clusters.R",
                      "remove_repeats.R",
                      "test_GRanges.R")
sapply(function_scripts, function(path){source(file = paste0("functions/", path))})


#Using a keep_cols list, remove unwanted metadata from GRanges
condense_metadata <- function(sites, keep_cols){
  tot_metadata <- names(mcols(sites))
  keep_positions <- match(keep_cols, tot_metadata)
  cleaned_sites <- sites[, keep_positions]
  return(cleaned_sites)
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

#intSiteStandardize, change name at some point
.intSiteStandardize <- function(sites, windowSize=5L, grouping=1, condense=FALSE, 
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


#Find all sites/clones using a list of position ID's (posid)
find_sites <- function(sites, posID){
  request <- do.call(c, lapply(1:length(posID), function(i){
    sites[sites$posID == posID[[i]],]
  }))
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