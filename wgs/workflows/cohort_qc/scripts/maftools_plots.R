#!/usr/bin/env Rscript

library(maftools)


oncoplot = function(read_maf, oncoplot_path, genes){
    png(filename=oncoplot_path, units="px", width=1600, height=1600, res=300)
    
    maftools::oncoplot(maf=read_maf, genes=genes)
    dev.off()
}

maf_summary = function(read_maf, maf_summary_path){
    png(filename=maf_summary_path, units="px", width=1600, height=1600, res=300)
    maftools::plotmafSummary(maf = read_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
}

somatic_interaction_plot = function(read_maf, somatic_interaction_plot_path){
    png(filename=somatic_interaction_plot_path, units="px", width=1600, height=1600, res=300)
    maftools::somaticInteractions(maf = read_maf, top=15, pvalue = c(0.05, 0.1))
    dev.off()
}


main = function(){
    args = commandArgs(trailingOnly=TRUE)
    genes=c("PPM1D", "TP53",  "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2","ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC")

    maf_file = args[1]
    cn=args[2]
    filtered_maf_file = args[3]
    oncoplot_path = args[4]
    somatic_interactions = args[5]

    mafsummary = args[6]
    vcNames = c("Frame_Shift_Del_somatic", "Frame_Shift_Ins_somatic", "Splice_Site_somatic", "Translation_Start_Site_somatic",
"Nonsense_Mutation_somatic", "Nonstop_Mutation_somatic", "In_Frame_Del_somatic","In_Frame_Ins_somatic", "Missense_Mutation_somatic",
"Frame_Shift_Del_germline", "Frame_Shift_Ins_germline", "Splice_Site_germline", "Translation_Start_Site_germline",
"Nonsense_Mutation_germline", "Nonstop_Mutation_germline", "In_Frame_Del_germline","In_Frame_Ins_germline", "Missense_Mutation_germline")
    unfiltmaf = maftools::read.maf(maf=maf_file, cnTable=cn, vc_nonSyn=vcNames)
    filtmaf = maftools::read.maf(maf=filtered_maf_file, cnTable=cn, vc_nonSyn=vcNames)
    
    oncoplot(filtmaf, oncoplot_path, genes)
    somatic_interaction_plot(unfiltmaf, somatic_interactions)
    maf_summary(unfiltmaf, mafsummary)

}

main()
