library(TitanCNA)
version <- "0.1.2"
args <- commandArgs(TRUE)

tumWig <- args[1]
normWig <- args[2]
gc <- args[3]
map <- args[4]
target_list <- args[5]
outfile <- args[6]
genometype <- args[7]

message('titan: Correcting GC content and mappability biases...')

if( target_list!="NULL" ){
    target_list_data <- read.table(target_list, sep="\t", header=F, stringsAsFactors=F)
    cnData <- correctReadDepth(tumWig, normWig, gc, map, genomeStyle=genometype, targetedSequence = target_list_data)
} else {
    cnData <- correctReadDepth(tumWig, normWig, gc, map, genomeStyle=genometype)
}

write.table(cnData, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")
