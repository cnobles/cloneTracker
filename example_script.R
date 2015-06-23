source("cloneTracker.utils.R")

data_names <- c("gr1", "gr2", "gr3", "gr4", "gr5")

load("ex_Data/gr1.RData")
load("ex_Data/gr2.RData")
load("ex_Data/gr3.RData")
load("ex_Data/gr4.RData")
load("ex_Data/gr5.RData")

gr3 <- intSiteCollapse(gr3)
gr4 <- intSiteCollapse(gr4)

hits_gr1xgr2 <- cloneTracker(query = gr1, subject = gr2, return_hits = 3)

grl <- GRangesList(gr1, gr2, gr3, gr4, gr5)
names(grl) <- data_names

grl <- lapply(grl, intSiteCollapse)

hits_list <- cloneTrackerTotal(grl, return_hits = 3)
names(hits_list) <- cloneTrackerTotalNames(site_names = data_names)

total_hits <- reduce(do.call(c, lapply(1:length(hits_list), 
                                       function(i){hits_list[[i]]})))

total_hits <- generate_posid(sites = total_hits)

logic_table <- generate_logic_table(hits = total_hits, list = grl)

