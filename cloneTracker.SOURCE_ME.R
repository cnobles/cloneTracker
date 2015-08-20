#Dependancies check, if required packages are not loaded, script is haulted
dependancies <- c("dplyr", "IRanges", "GenomicRanges", "igraph", "devtools")

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
function_scripts <- c("/track_clones.R", 
                      "/determine_abundance.R", 
                      "/standardize_intsites.R", 
                      "/normalize_multihit_clusters.R",
                      "/remove_repeats.R",
                      "/db_to_granges.R",
                      "/condense_intsites.R",
                      "/cloneTracker.utils.R",
                      "/serial_cluster.R",
                      "/normalize_intsite_positions.R",
                      "/ex_Data/test_GRanges.R")
raw.script.url <- "https://raw.githubusercontent.com/cnobles/cloneTracker/master"
sapply(function_scripts, function(path){source_url(paste0(raw.script.url, path))})
message("cloneTracker functions loaded.")
rm(raw.script.url, function_scripts)
