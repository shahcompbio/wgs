library(maftools)


oncoplot = function(read_maf, oncoplot_path, genes){
    png(filename=oncoplot_path, units="px", width=1600, height=1600, res=300)
    
    maftools::oncoplot(maf=read_maf, genes=genes)
    dev.off()
}

maf_summary = function(read_maf, maf_summary_path){
        png(filename=maf_summary_path, units="px", width=1600, height=1600, res=300)
        plot.new()
        dev.off()
    try({ 
        png(filename=maf_summary_path, units="px", width=1600, height=1600, res=300)
        maftools::plotmafSummary(maf = read_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
        dev.off()
    })
}

somatic_interaction_plot = function(read_maf, somatic_interaction_plot_path){
        png(filename=somatic_interaction_plot_path, units="px", width=1600, height=1600, res=300)
        plot.new()
        dev.off()
    try({ 
        png(filename=somatic_interaction_plot_path, units="px", width=1600, height=1600, res=300)
        maftools::somaticInteractions(maf = read_maf, top=15, pvalue = c(0.05, 0.1))
        dev.off()
    })
}


main = function(){
    args = commandArgs(trailingOnly=TRUE)
    genes=c("PPM1D", "TP53",  "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2","ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC")

    maf_file = args[1]
    vcNames=args[2]
    cn=args[3]
    oncoplot_path = args[4]
    somatic_interactions = args[5]

    mafsummary = args[6]

    vcNames=read.table(vcNames,header=TRUE)$Variant_Classification

    maf = maftools::read.maf(maf=maf_file, cnTable=cn, vc_nonSyn=vcNames)
    
    oncoplot(maf, oncoplot_path, genes)
    maf_summary(maf, mafsummary)
    somatic_interaction_plot(maf, somatic_interactions)


}


main()
