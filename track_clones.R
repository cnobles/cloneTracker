#Return all clones present in more than 1 set of data from a list of sites

track_clones <- function(sites.list, gap=5L, track.origin=TRUE){
  grl.sites <- sites.list
  
  if(class(sites.list) == "list"){grl.sites <- GRangesList(sites.list)}
  
  if(track.origin){
    grl.sites <- GRangesList(lapply(1:length(sites.list), function(i){
      sites <- sites.list[[i]]
      sites$origin <- rep(names(sites.list[i]), length(sites))
      sites
    }))
  }
  
  ovlp.grps <- findOverlaps(grl.sites, maxgap = gap, ignoreSelf = TRUE, ignoreRedundant = TRUE)
  if(length(ovlp.grps) > 0){
   ovlp.sites <- unlist(GRangesList(lapply(1:length(ovlp.grps), function(i){
     query <- grl.sites[[queryHits(ovlp.grps[i])]]
     subject <- grl.sites[[subjectHits(ovlp.grps[i])]]
     hits <- findOverlaps(
       flank(query, -1, start = TRUE), 
       flank(subject, -1, start = TRUE), 
       maxgap = gap)
     if(length(hits) > 0){
       sites <- c(query[queryHits(hits)], subject[subjectHits(hits)])
     }else{
       sites <- GRanges()
     }
     sites  
   })))
  }else{
   message("No overlaping sites found between groups.")
  }
  
  if(length(ovlp.grps) > 0){
    ovlp.sites$posid <- generate_posID(ovlp.sites)
    sites.dfr <- distinct(as.data.frame(ovlp.sites, row.names = NULL))
    ranges <- IRanges(start = sites.dfr$start, end = sties.dfr$end)
    sites.gr <- GRanges(
      seqnames = sites.dfr$seqnames,
      ranges = ranges,
      strand = sites.dfr$stand
    )
    mcols(sites.gr) <- sites.dfr[,c(6:length(sites.dfr))]
    ovlp.list <- split(sites.gr, sites.gr$posid)
  }else{
    ovlp.list <- GRangesList()
  }
  message(paste("Number of overlaping sites:", length(ovlp.list)))
  ovlp.list
}
