#Return all clones present in more than 1 set of data from a list of sites

track_clones <- function(sites.list, gap=5L, track.origin=TRUE, standardize=TRUE, ...){
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
     hits <- findOverlaps(query, subject, maxgap = gap)
     c(query[queryHits(hits)], subject[subjectHits(hits)])
   })))
  }else{
   message("No overlaping sites found between groups.")
   std.list <- GRangesList()
  }
  
  if(length(ovlp.grps) > 0){
   if(standardize){
     std.sites <- standardize_intsites(ovlp.sites, std.gap = 1L, standardize_breakpoints = FALSE)
   }else{
     std.sites <- ovlp.sites
   }
  
   std.sites$posid <- generate_posID(std.sites)
   std.list <- split(std.sites, std.sites$posid)
  }
  
  std.list
}
