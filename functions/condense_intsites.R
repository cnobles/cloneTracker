condense_intsites <- function(sites_to_condense, grouping = NULL, return.abundance = FALSE, 
                              method = "fragLen", replicates = "replicates"){
  mcols <- mcols(sites_to_condense)
  
  grp <- which(grouping == names(mcols))
  if(length(grp) != 0){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites_to_condense))
  }
  
  sites_to_condense$group <- group
  
  sites.list <- split(sites_to_condense, sites_to_condense$group)
  
  condensed.grl <- GRangesList(lapply(sites.list, function(sites){
    sites <- flank(sites, -1, start = TRUE)
    sites$posID <- generate_posID(sites) 
    first.hit <- match(unique(sites$posID), sites$posID)
    condensed.sites <- sites[first.hit]
    condensed.sites
  }))

  condensed.gr <- unlist(condensed.grl)
  
  if(return.abundance){
    abund.dfr <- determine_abundance(sites_to_condense, grouping = grouping, 
                                     replicates = replicates, method = method)
    mcols <- as.data.frame(mcols(condensed.gr))
    mcols.list <- split(mcols, mcols$group)
    abund.list <- split(abund.dfr, abund.dfr$group)
    mcols.list <- mcols.list[names(abund.list)]
    all.cols <- bind_rows(lapply(1:length(mcols.list), function(i){
      merge(mcols.list[[i]], abund.list[[i]], by="posID")
    }))
    order <- match(generate_posID(condensed.gr), all.cols$posID)
    all.cols <- all.cols[order,]
    mcols(condensed.gr) <- all.cols
    condensed.gr$group.x <- NULL
    condensed.gr$group.y <- NULL
  }
  names(condensed.gr) <- NULL
  message("Check to make sure all metadata columns are still relavent.")
  condensed.gr
}