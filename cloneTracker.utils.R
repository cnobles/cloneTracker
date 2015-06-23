#Dependancies check, if required packages are not loaded, script is haulted
dependancies <- c("IRanges", "GenomicRanges")

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


#Utility functions for clone tracking using integration sites as markers
#Remove width data from called intsites, returned data is used for
#clone marking
intSiteCollapse <- function(sites, use_names=TRUE){
  strand <- as.vector(strand(sites))
  adj_start <- ifelse(strand == "+", start(sites), end(sites))
  adj_width <- rep(1, length(sites))
  names <- ifelse(use_names, names(sites), NULL)
  adj_ranges <- IRanges(start = adj_start, width = adj_width, names = names)
  adj_sites <- GRanges(seqnames = seqnames(sites), 
                       ranges = adj_ranges, strand = strand(sites))
  mcols(adj_sites) <- mcols(sites)
  return(adj_sites)
}

#Return clones present in query and subject with default conditions
cloneTracker <- function(query, subject, maxgap=5L, minoverlap=1L, 
                            type="any", select="all", return_hits=1, ...){
  
  overlap_check_table <- table(
    countOverlaps(query = query, subject = subject, maxgap = maxgap, 
                  minoverlap = minoverlap, type = type)
    )
  
  if(overlap_check_table[[1]] == length(query)){
    empty_grange <- GRanges()
    return(empty_grange)
  }
  
  overlaps <- findOverlaps(query = query, subject = subject, 
                           maxgap = maxgap, minoverlap = minoverlap, 
                           type = type, select = select)
  
  query_hits <- do.call(c, sapply(1:length(overlaps), function(i){
    hit <- query[queryHits(overlaps)[i]]
  }))
  names(query_hits) <- 1:length(query_hits)
  subject_hits <- do.call(c, sapply(1:length(overlaps), function(i){
    hit <- subject[subjectHits(overlaps)[i]]
  }))
  names(subject_hits) <- 1:length(subject_hits)
  total_hits <- c(query_hits, subject_hits)
  names(total_hits) <- 1:length(total_hits)
  hits_list <- list(query_hits, subject_hits, total_hits)
  return(hits_list[[return_hits]])
}

#Return all clones present in more than 1 set of data from a list of sites
#Function makes every pairwise comparison posible from the list given
cloneTrackerTotal <- function(sites, maxgap=5L, minoverlap=1L, 
                                  type="any", select="all", 
                                  return_hits=3, ...){
  compar_matrix <- combn(1:length(sites), 2)
  total_hits <- lapply(1:ncol(compar_matrix), function(i){
    query <- sites[[ compar_matrix[[1,i]] ]]
    subject <- sites[[ compar_matrix[[2,i]] ]]
    
    overlap_check_table <- table(
      countOverlaps(query = query, subject = subject, 
                    maxgap = maxgap, minoverlap = minoverlap, 
                    type = type)
      )
    
    if(overlap_check_table[[1]] == length(query)){
      empty_grange <- GRanges()
      return(empty_grange)
    }
    
    overlaps <- findOverlaps(query = query, subject = subject, 
                             maxgap = maxgap, minoverlap = minoverlap, 
                             type = type, select = select)
    
    query_hits <- do.call(c, sapply(1:length(overlaps), function(j){
      hit <- query[queryHits(overlaps)[j]]
    }))
    names(query_hits) <- 1:length(query_hits)
    
    subject_hits <- do.call(c, sapply(1:length(overlaps), function(j){
      hit <- subject[subjectHits(overlaps)[j]]
    }))
    names(subject_hits) <- 1:length(subject_hits)
    
    total_hits <- c(query_hits, subject_hits)
    names(total_hits) <- 1:length(total_hits)
    hits_list <- list(query_hits, subject_hits, total_hits)
    return(hits_list[[return_hits]])
  })
  return(total_hits)
}

#Function generates the names for the overlaping comarisons for 
#return_total_overlaps
cloneTrackerTotalNames <- function(site_names){
  compar_matrix <- combn(site_names, 2)
  total_names <- sapply(1:ncol(compar_matrix), function(i){
    name <- paste0(compar_matrix[[1,i]], "_X_", compar_matrix[[2,i]])
  })
  return(total_names)
}

#Find all sites/clones using a list of position ID's (posid)
find_sites <- function(sites, posid){
  request <- do.call(c, lapply(1:length(posid), function(i){
    sites[sites$posid == posid[[i]],]
  }))
}

#After clustering sites by hiReadsProcessor::clusterSites, select only
#top hits, giving the dominant site for the cluster
collect_top_hits_for_int_sites <- function(sites){
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
  sites$posid <- paste0(chr, strand, pos)
  return(sites)
}

#A logic table simply reports on the presence or absence of a clone in a 
#data set
generate_logic_table <- function(hits, list){
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
generate_abund_table <- function()
  
  
#Though likely modified for specific situations, score based on logic table
#to determine best ranking or tracked clones
score_tracking_hits <- function(i, table){
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
  
  return(score)
}