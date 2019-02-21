library(IchorCNA)
version <- "0.1.2"
args <- commandArgs(TRUE)

tumWig <- args[1]
normWig <- args[2]
gc <- args[3]
map <- args[4]
target_list <- args[5]
outfile <- args[6]
genometype <- args[7]

