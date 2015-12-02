#Test GRange
test <- GRanges(seqnames = Rle(values = c("chr3"),
                               lengths = c(15)),
                ranges = IRanges(start = c(rep(6,10),rep(7,3),8,10),
                                 end = c(25,20,25,25,24,20,19,30,30,25,24,24,20,29,23)),
                strand = Rle(values = c("+"),
                             lengths = c(15)))
