#Population calculation functions

#Calculate the Shannon Diversity index given a set of frequencies / proportions
calc_shannon <- function(x, base = exp(1)){
  x <- x/sum(x, na.rm = TRUE)
  shannon <- -x*log(x, base)
  shannon <- sum(shannon, na.rm = TRUE)
  shannon
}

#Calculate Gini Index for a sample given frequencies / proportions
calc_gini <- function(x){
  stopifnot(require(reldist))
  x <- x/sum(x, na.rm = TRUE)
  gini <- gini(x)
  gini
}

#Calculate the Entropy of a sample set as defined by Adaptive Biotechnologies
calc_entropy <- function(x){
  x <- x/sum(x, na.rm = TRUE)
  entropy <- -sum(x*log(x, base = 2))
  entropy
}

#Calculate the clonality as defined by Adaptive Biotechnologies for TCR and IgH sequence analysis
calc_clonality <- function(x){
  x <- x/sum(x, na.rm = TRUE)
  clonality <- 1+sum(x*log(x, base = length(x)))
  clonality
}