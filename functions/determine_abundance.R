#requires dplyr::distinct and sonicLength::estAbund
intSiteAbund <- function(sites, method="fragLen"){
  sites$posID <- generate_posID(sites)
  mcols <- mcols(sites)
  mcols(sites) <- NULL
  
  reps <- which("replicates" == names(mcols))
  if(length(reps) != 0){
    replicates <- mcols$replicates
  }else{
    replicates <- rep("group1", length(sites))
  }
  
  posID <- generate_posID(sites)
  fragLen <- width(sites)
  sites.dfr <- data.frame("posID"=posID, 
                          "fragLen"=fragLen, 
                          "replicates"=replicates)
  
  if(method == "fragLen"){
    abundCalc <- function(locations, fragLen, replicates){
      if(length(unique(sites.dfr$replicates)) == 1){
        locationID <- locations
      }else{
        locationID <- paste0(replicates, ":", locations)
      }
      
      dfr <- data.frame("locationID"=locationID, "fragLen" = fragLen)
      dfr_dist <- distinct(dfr)
      sites_list <- split(dfr_dist, dfr_dist$locationID)
      abundances <- sapply(sites_list, function(x){nrow(x)})
      abundances <- abundances[unique(dfr$locationID)]
      
      if(length(unique(sites.dfr$replicates)) == 1){
        posID <- names(abundances)
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances)
      }else{
        group_posID <- strsplit(names(abundances), split = ":")
        group <- sapply(1:length(group_posID), function(i){group_posID[[i]][1]})
        posID <- sapply(1:length(group_posID), function(i){group_posID[[i]][2]})
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances,
                                "replicates" = group)
      }
      abund.dfr
    }
  }else if(method == "estAbund"){
    abundCalc <- function(locations, fragLen, replicates=NULL){
      if(length(unique(sites.dfr$replicates)) == 1){
        theta_list <- estAbund(locations=locations, lengths=fragLen)
      }else{
        theta_list <- estAbund(locations=locations, lengths=fragLen, replicates=replicates)
      }
      posID <- names(theta_list$theta)
      abundances <- theta_list$theta
      if(length(unique(sites.dfr$replicates)) == 1){
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances)
      }else{
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances,
                                "replicates" = theta_list$data$replicates)
      }
      abund.dfr
    }
  }else{
    stop("Must choose either fragLen or estAbund for method.")
  }
  
  abund.dfr <- abundCalc(locations = sites.dfr$posID,
                         fragLen = sites.dfr$fragLen, 
                         replicates = sites.dfr$replicates)
  
  abund.dfr$estAbund <- round(abund.dfr$estAbund)
  abund.dfr$estAbundProp <- abund.dfr$estAbund/sum(abund.dfr$estAbund)
  abund.dfr$estAbundRank <- rank(-1*abund.dfr$estAbundProp, ties.method="max")
  rownames(abund.dfr) <- NULL
  abund.dfr
}
