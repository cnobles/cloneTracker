#Test GRange
test <- GRanges(seqnames = Rle(values = c("chr1", "chr2", "chr3"),
                               lengths = c(5, 5, 5)),
                ranges = IRanges(start = c(2,7,10,7,8,15,15,17,16,14,6,7,6,5,2),
                                 end = c(18,25,25,28,24,30,30,30,29,30,16,20,19,16,25)),
                strand = Rle(values = c("+","-","+"),
                             lengths = c(5,5,5)))
