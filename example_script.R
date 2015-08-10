source("cloneTracker.SOURCE_ME.R")

data_names <- c("gr1", "gr2", "gr3", "gr4", "gr5")

load("ex_Data/gr1.RData")
load("ex_Data/gr2.RData")
load("ex_Data/gr3.RData")
load("ex_Data/gr4.RData")
load("ex_Data/gr5.RData")

grl <- GRangesList(gr1, gr2, gr3, gr4, gr5)
names(grl) <- data_names

grl <- lapply(grl, function(x){flank(x, -1, start = TRUE)})

cloneList <- track_clones(grl)
