condense_intsites <- function(sites_to_condense, grouping = NULL, return.abundance = FALSE, 
                              method = "fragLen", replicates = "replicates"){
  mcols <- mcols(sites_to_condense)
  
  grp <- which(grouping == names(mcols))
  if(length(grp) != 0){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites_to_condense))
  }
  
  sites_to_condense$posID <- generate_posID(sites_to_condense)
  sites_to_condense$groupID <- paste0(group, "^", sites_to_condense$posID)
  
  groupIDs <- unique(sites_to_condense$groupID)
  first.hits <- match(groupIDs, sites_to_condense$groupID)
  
  condensed.gr <- flank(sites_to_condense, -1, start = TRUE)
  condensed.gr <- condensed.gr[first.hits]
  
  if(return.abundance){
    abund.dfr <- determine_abundance(sites_to_condense, grouping = grouping, 
                                     replicates = replicates, method = method)
    abund.dfr$groupID <- paste0(abund.dfr$group, "^", abund.dfr$posID)
    mcols <- mcols(condensed.gr)
    all.cols <- merge(mcols, abund.dfr, by = "groupID")
    order <- match(condensed.gr$groupID, all.cols$groupID)
    all.cols <- all.cols[order,]
    mcols(condensed.gr) <- all.cols
    condensed.gr$posID.x <- NULL
    condensed.gr$posID.y <- NULL
    condensed.gr$posID <- generate_posID(condensed.gr)
  }
  condensed.gr$groupID <- NULL
  condensed.gr$group <- NULL
  names(condensed.gr) <- NULL
  message("Check to make sure all metadata columns are still relavent.")
  condensed.gr
}