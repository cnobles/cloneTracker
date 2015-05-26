source("cloneTracker.utils.R")

load("ex_Data/gr1.RData")
load("ex_Data/gr2.RData")

gr1 <- collapse_ranges_to_intSites(gr1)
gr2 <- collapse_ranges_to_intSites(gr2)

hits_gr1xgr2 <- return_overlaps(query = gr1, subject = gr2, return_hits = 3)


load("ex_Data/gr3.RData")
load("ex_Data/gr4.RData")

gr_list <- list(gr1, gr2, gr3, gr4)
names_gr_list <- c("gr1", "gr2", "gr3", "gr4")
names(gr_list) <- names_gr_list

gr_list <- lapply(gr_list, collapse_ranges_to_intSites)

hits_list <- return_total_overlaps(gr_list, return_hits = 3)
names(hits_list) <- names_total_overlaps(site_names = names_gr_list)

total_hits <- reduce(do.call(c, lapply(1:length(hits_list), 
                                       function(i){hits_list[[i]]})))

total_hits <- generate_posid(sites = total_hits)

logic_table <- generate_logic_table(hits = total_hits, list = gr_list)

