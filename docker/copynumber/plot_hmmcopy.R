#!/usr/bin/env Rscript
library(HMMcopy)

args <- commandArgs(TRUE)

tumour_copy <- args[1]
hmmcopy_res <- args[2]
correction_plots_dir <- args[3]
bias_plots <- args[4]
hmmcopy_plots_dir <- args[5]


dir.create(hmmcopy_plots_dir, recursive=TRUE)
dir.create(correction_plots_dir, recursive=TRUE)
dir.create(dirname(bias_plots), recursive=TRUE)

load(tumour_copy)
tumour_copy <- infile_copy

load(hmmcopy_res)
hmmcopy_res <- tumour_segments

bias_fig <- bias_plots
bias_fig <- sub(".csv", ".pdf", bias_plots)

#plot the effects of GC/mappability correction

pdf(bias_fig)
par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
plotBias(tumour_copy, pch = 20, cex = 0.5)
dev.off()



chr<-unique(space(tumour_copy))
for (i in chr){
    correction_fig <- paste(correction_plots_dir,"/chr_",i,".jpg",sep="")

    jpeg(correction_fig, units="in", width=11, height=8.5, res=500)
    par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))

    tryCatch({
        plotCorrection(tumour_copy, chr=i)
    },
    error=function(err){
    print(paste('Warning: Correction plot skipped for chromsome',i,sep = " "))
    #print(err)
    })

    dev.off()
}


# plot copy numbers
chr <- unique(space(infile_copy))

for (i in chr){
    outplot <- paste(hmmcopy_plots_dir,"/chr_",i,".pdf",sep="")
    pdf(outplot)
    par(mfrow = c(1, 1))
    par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 2, 2, 2), mgp = c(1, 0.5, 0))

    plotSegments(tumour_copy, hmmcopy_res, chr = i, pch = ".", ylab = "Copy Number", xlab = "Chromosome Position", main=paste("Chr ",i,sep=""))
    cols <- stateCols()
    legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"), fill = cols, horiz = TRUE, bty = "n", cex = 0.5)
    dev.off()
}



