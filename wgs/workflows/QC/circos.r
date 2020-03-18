library(circlize)
library(ComplexHeatmap)

#add a copy number track to the circos plot
#copy_number -> copy number data.frame
#gene_annotations -> gene_annotations data.frame
plot_copy_number = function(copy_number, gene_annotations){
  require(circlize)

  cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
  names(cnCol) <- c(0:500)


  circos.track(factors= paste0("chr", copy_number$Chr), ylim=c(-4, 6), x=copy_number$Position, "track.height" = 0.15,
               panel.fun = function(x, y) {
                 for (state in -4:6){
                  circos.lines(c(0, max(x)), c(state, state), col = "#00000040")
                 }
               })

  for (state in unique(copy_number$state)){
    x = copy_number$Position[copy_number$state == state]
    y = copy_number$LogRatio[copy_number$state == state]
    f = copy_number$Chr[copy_number$state == state]
    col = cnCol[state+1]
    circos.trackPoints(factors= paste0("chr", f),  x = x, y = y, cex = 0.25, col=col)
  }

  apply(gene_annotations,1,function(row)
        circos.trackLines(factors=rep(paste0("chr",strtoi(row["chrom"])), 2),
                          x=rep(strtoi(row["pos"]), 2),
                          y=c(-4, 6),
                          col = row["color"])
        )
}

#add the svs to the circos plot as links
#svs -> data.frame of breakpoints
plot_svs = function(svs){
  require(circlize)
  for (sv in 1:nrow(svs)){
    chr1 = paste0("chr", svs[sv, "chromosome_1"])
    chr2 = paste0("chr", svs[sv, "chromosome_2"])
    p1 =  svs[sv, "position_1"]
    p2 =  svs[sv, "position_2"]
    color = svs[sv, "color"]
    circos.link(chr1, p1, chr2, p2, col = color, h.ratio=0.7)


  }
}

#add a legend to the plot
add_legend = function(){
  cn_leg = Legend(labels = c("LogRatio"), type = "points",
                  legend_gp =c("black"), title_position = "topleft", title = "Copy Number")

  svs = Legend(labels = c('foldback', 'unbalanced', 'duplication','deletion', 'inversion', 'balanced'), type = "points",
               legend_gp = gpar(col =c(2, 3, 4, 1, 6, 8), lwd = 2), title_position = "topleft",
               title = "Structural Variants")

  labels = c("CCNE1","ERBB2",  "KRAS",   "MYC", "PIK3CA", "MECOM",  "RB1",
             "PTEN",   "BRCA1",  "BRCA2",  "RAD51C", "PALB2")

  cols = c("blue", "green", "red", "cyan", "magenta",  "yellow", "lightsteelblue", "tan","black",
           "grey","saddlebrown", "pink")
  gene_annotations = Legend(labels = labels, type = "points",
                            legend_gp = gpar(col = cols, lwd = 2), title_position = "topleft", title = "Gene Annotations")

  leg = packLegend(cn_leg, svs, gene_annotations)

  draw(leg,x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
  #draw(gene_annotations,x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "top"))


}

#initialize a circos plot and fill with copy number ring, svs and legend
#copy_number -> copy number data.frame
#annotations -> gene annotations data.frame
#svs -> sv data.frame
circos_plot = function(copy_number, annotations,svs){
  require(circlize)



  circos.initializeWithIdeogram(
    chromosome.index = paste0("chr", c(1:22, "X")),
  )

  plot_copy_number(copy_number, annotations)
  plot_svs(svs)
  circos.clear()

  add_legend()
}

#make a pdf and write circos plot to it
#copy_number -> copy number data.frame
#annotations -> gene annotations data.frame
#svs -> sv data.frame
#outfile -> path for output pdf file
make_circos = function(copy_number, annotations, svs, outfile){

  #save plot to pdf
  pdf(outfile, 8, 13)

  #create plot
  circos_plot(copy_number, annotations, svs)

  #turn off pdf
  dev.off()

}


#example inputs
#example inputs
anns = "/Users/abramsd/work/CODE/QC-shipped/SHIPPED/CIRCOS/annotations.tsv"
anns = read.csv(anns, sep="\t")

cn7 = "/Users/abramsd/work/CODE/QC-shipped/SHIPPED/007_binned_cn.csv"
cn7 = read.csv(cn7, sep="\t")
svs7 = "/Users/abramsd/work/CODE/QC-shipped/SHIPPED/007svs"
svs7 = read.csv(svs7, sep="\t")

make_circos(cn7, anns, svs7, "example_output.pdf")
