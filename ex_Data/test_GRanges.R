#Test GRange
test <- GRanges(seqnames = Rle(values = c("chr1", "chr2", "chr3"),
                               lengths = c(5, 5, 5)),
                ranges = IRanges(start = c(2,7,10,7,8,15,10,5,16,17,6,8,6,5,4),
                                 end = c(18,25,25,28,24,30,30,30,29,28,16,20,17,16,20)),
                strand = Rle(values = c("+","-","+"),
                             lengths = c(5,5,5)))
