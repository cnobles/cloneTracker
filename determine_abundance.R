#requires dplyr::distinct and sonicLength::estAbund
determine_abundance <- function(sites, grouping=NULL, replicates=NULL, method="fragLen"){
  sites$posID <- generate_posID(sites)
  mcols <- mcols(sites)
  mcols(sites) <- NULL
  
  grp <- which(grouping == names(mcols))
  if(length(grp) != 0){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites))
  }
  
  reps <- which(replicates == names(mcols))
  if(length(reps) != 0){
    replicates <- as.character(mcols[,reps])
  }else{
    replicates <- rep("1", length(sites))
  }
  
  sites.dfr <- data.frame("posID"= generate_posID(sites), 
                          "fragLen"= width(sites), 
                          "replicates"=replicates,
                          "group"=group)
  
  if(method == "fragLen"){
    abundCalc <- function(locations, fragLen, replicates, group){
      group <- unique(group)
      dfr <- data.frame("locationID" = locations, 
                        "fragLen" = fragLen,
                        "replicates" = replicates)
      dfr_dist <- distinct(dfr)
      sites_list <- split(dfr_dist, dfr_dist$locationID)
      abundances <- sapply(sites_list, function(x){nrow(x)})
      abundances <- abundances[unique(dfr$locationID)]
      
      abund.dfr <- data.frame("posID" = names(abundances),
                              "estAbund" = abundances,
                              "group" = rep(group, length(abundances)))
      abund.dfr
    }
  }else if(method == "estAbund"){
    abundCalc <- function(locations, fragLen, replicates, group){
      group <- unique(group)
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
                                "replicates" = length(unique(theta_list$data$replicates)))
      }
      abund.dfr$group <- group
      abund.dfr
    }
  }else{
    stop("Must choose either fragLen or estAbund for method.")
  }
  
  sites.list <- split(sites.dfr, sites.dfr$group)
 
  abund.list <- lapply(sites.list, function(x){
    abundCalc(locations = x$posID,
              fragLen = x$fragLen, 
              replicates = x$replicates,
              group = x$group)
  })
  
  abund.list <- lapply(abund.list, function(x){
    x <- x %>%
      mutate(totalAbund = sum(x$estAbund)) %>%
      mutate(estAbund = round(estAbund)) %>%
      mutate(relAbund = estAbund/totalAbund) %>%
      mutate(relRank = rank(-1*relAbund, ties.method="max"))
    x
  })
  
  abund.dfr <- bind_rows(abund.list)
  
  abund.dfr$totalAbund <- NULL
    
  abund.dfr <- abund.dfr[, c("posID", "group", "estAbund", "relAbund", "relRank")]
  abund.dfr
}
