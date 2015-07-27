#Test GRange
test <- GRanges(seqnames = Rle(values = c("chr1", "chr2", "chr3"),
                               lengths = c(5, 10, 7)),
                ranges = IRanges(start = c(2,7,10,7,8,15,15,17,16,14,14,6,7,6,5,2,5,8,11,14,17,20),
                                 width = rep(30, 22)),
                strand = Rle(values = c("+","-","+"),
                             lengths = c(8,7,7)))
