#!/usr/bin/env Rscript
library(HMMcopy)

args <- commandArgs(TRUE)

infile <- args[1]
gc <- args[2]
map <- args[3]
mapcutoff <- as.numeric(args[4])
outfile <- args[5]
outobj <- args[6]

# correct read depths
infile_copy <- correctReadcount(wigsToRangedData(infile, gc,map), mappability=mapcutoff)

save( infile_copy , file=paste(outobj))


infile_copy <- as.data.frame(infile_copy)
colnames(infile_copy) <- c("chr","start","end","width","reads","gc","map","valid","ideal","cor.gc","cor.map","copy")
write.table(infile_copy, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")
