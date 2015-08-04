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
  
  condensed.list <- lapply(sites.list, function(sites){
    sites <- flank(sites, -1, start = TRUE)
    sites$posID <- generate_posID(sites) 
    first.hit <- match(unique(sites$posID), sites$posID)
    condensed.sites <- sites[first.hit]
    condensed.sites
  })
  #  #Make sure sites have been standardized
  #  sites.reduced <- flank(sites, -1, start=TRUE)
  #  sites.reduced <- unlist(reduce(sites.reduced, min.gapwidth=0L, 
  #                                 with.revmap=TRUE))
  #  sites.reduced$posID <- generate_posID(sites.reduced)
  #  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
    #Order original sites by revmap  
  #  condensed.sites <- sites[unlist(sites.reduced$revmap)]
  #  condensed.sites$posID <- generate_posID(condensed.sites)
  
  #  if(length(sites) > 0){
  #    condensed.sites <- split(condensed.sites, 
  #                             Rle(values = seq(length(sites.reduced)), 
  #                                 lengths = sites.reduced$counts))
  #  }
  
  #  #Condense the reads with same standardized starts and remove breakpoint info
  #  #Wrapped by the collection of metadata columns
  #  metacols <- do.call(rbind, lapply(condensed.sites, function(x){
  #    mcols(x[1])
  #  }))
  #  condensed.sites <- unlist(reduce(flank(condensed.sites, -1, start=TRUE)))
  #  mcols(condensed.sites) <- merge(mcols(sites.reduced), metacols, by="posID")
  #  condensed.sites
  #})
  
  condensed.gr <- do.call(c, lapply(1:length(condensed.list), function(i){condensed.list[[i]]}))
  
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
    mcols(condensed.gr) <- all.cols
    condensed.gr$group.x <- NULL
    condensed.gr$group.y <- NULL
  }
  condensed.gr$revmap <- NULL
  condensed.gr$counts <- NULL
  condensed.gr$group <- NULL
  message("Check to make sure all metadata columns are still relavent.")
  condensed.gr
}