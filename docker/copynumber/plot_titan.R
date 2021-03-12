#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

version <- "0.1.4"

library(TitanCNA)
library(SNPchip)  ## use this library to plot chromosome idiogram

obj_file <- args[1]
outdir <- args[2]
numClusters <- args[3]
chr <- args[4]
chr <- eval(parse(text=chr))
ploidy <- args[5]

# load a workspace into the current session
load(obj_file)

if(ploidy == "NULL"){
    ploidy = convergeParams$ploidy
}

if (numClusters == "NULL"){
    numClusters = convergeParams$numClusters
}

#### PLOT RESULTS ####
finoutdir <- paste0(outdir, "/cluster_", numClusters, "_ploidy_", ploidy, "/")

dir.create(finoutdir, recursive=TRUE)

norm <- tail(convergeParams$n, 1)
ploidy <- tail(convergeParams$phi, 1)

for (i in chr) {
    outplot <- paste0(finoutdir, "cluster_", numClusters, "_chr", i, ".png")
    png(outplot, width=1200, height=1000, res=100, type="cairo")

    # if 2 or fewer clusters, use c(4,1) panels
    if (as.numeric(numClusters) <= 2) {
		par(mfrow=c(4,1))
	} else {
		par(mfrow=c(3,1))
	}

    ## PLOT LOG RATIO (CNA) ##
    tryCatch({
        plotCNlogRByChr(results, i, ploidy=ploidy, geneAnnot=NULL, spacing=4,
                        norm=norm,ylim=c(-4,6), cex=0.5, main=paste("Chr ",i,sep=""), xlab="")
    },
    error=function(err){
        print(paste('Warning: plotCNlogRByChr skipped for chromosome', i, 'due to', err))
    })

    ## PLOT ALLELIC RATIOS (LOH) ##
    tryCatch({
        plotAllelicRatio(results, i, geneAnnot=NULL, spacing=4,
        				ylim=c(0,1), cex=0.5, main=paste("Chr ",i,sep=""), xlab="")
    },
    error=function(err){
        print(paste('Warning: plotAllelicRatio skipped for chromosome', i, 'due to', err))
    })

    ## PLOT CELLULAR PREVALENCE AND CLONAL CLUSTERS ##
    tryCatch({
        plotClonalFrequency(results, i, normal=norm, geneAnnot=NULL, spacing=4,
                        ylim=c(0,1), cex=0.5, main=paste("Chr ",i,sep=""), xlab="")
    },
    error=function(err){
        print(paste('Warning: plotClonalFrequency skipped for chromosome', i, 'due to', err))
    })

    if (as.numeric(numClusters) <= 2) {
        tryCatch({
            plotSubcloneProfiles(results, i, cex = 2, spacing=6, main=paste("Chr ",i,sep=""))
        },
        error=function(err){
            print(paste('Warning: plotSubcloneProfiles skipped for chromosome', i, 'due to', err))
        })
    }

    ## PLOT SUBCLONE PROFILE FOR 2 OR FEWER CLONAL CLUSTER RUNS ##
    tryCatch({
        if (as.numeric(numClusters) <= 2){
    		pI <- plotIdiogram(i, build="hg19", unit="bp", label.y=-4.25, new=FALSE, ylim=c(-2,-1))
    	} else {
    		pI <- plotIdiogram(i, build="hg19", unit="bp", label.y=-0.35, new=FALSE, ylim=c(-0.2,-0.1))
    	}

    },
    error=function(err){
        print(paste('Warning: plotIdiogram skipped for chromosome ', i, 'due to ', err))
    })

    dev.off()
}
