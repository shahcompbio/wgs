library(circlize)
library(ComplexHeatmap)

fill_cn_data_chroms = function(data, chrom_colname, expected_chroms){
  data_chroms = unique(data$chrom_colname)
  chroms_to_add = expected_chroms[! expected_chroms %in% data_chroms]
  for (chrom in chroms_to_add){
    r = nrow(data)+1
    data[r,"LogRatio"] = 0
    data[r,"Position"] = 0
    data[r,"state"] = 0
    data[r,"Chr"] = chrom

  }
  return(data)
}


#add a copy number track to the circos plot
#copy_number -> copy number data.frame
#gene_annotations -> gene_annotations data.frame
plot_copy_number = function(copy_number, gene_annotations){
  cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
  names(cnCol) <- c(0:500)
  chroms_needed = c(1:22, "X")
  if (! all (chroms_needed %in% copy_number$Chr)){
    copy_number = fill_cn_data_chroms(copy_number, "Chr", chroms_needed)
  }
  print( typeof(copy_number$Position))
  print(copy_number$Position)
  circos.track(factors= paste0("chr", copy_number$Chr), ylim=c(-2, 2), x=copy_number$Position, "track.height" = 0.15,
               panel.fun = function(x, y) {
                 for (state in -2:2){
                   circos.lines(c(0, max(x)), c(state, state), col = "#00000040")
                 }
               })

    for (state in unique(copy_number$state)){
    x = copy_number$Position[copy_number$state == state]
    y = copy_number$LogRatio[copy_number$state == state]
    f = copy_number$Chr[copy_number$state == state]
    col = cnCol[state+1]
    circos.trackPoints(factors= paste0("chr", f),  x = x, y = y, cex = 0.1, col=col)
  }

}

#add the svs to the circos plot as links
#svs -> data.frame of breakpoints
plot_svs = function(svs){
  if (! nrow(svs)){
    return()
  }
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
  col_fun = colorRamp2(c(1, 2, 3, 4, 5, 6, 7, 8, 9), c("#00FF00", "#006400", "#0000FF", "#880000",
                                                       "#BB0000", "#CC0000", "#DD0000", "#EE0000", "#FF0000"))

  cn_leg = Legend(labels = c("9+", "8",  "7",  "6",  "5",  "4",  "3",  "2",  "1"), type="points",
                  legend_gp=gpar(col=c("#FF0000", "#EE0000", "#DD0000", "#CC0000", "#BB0000", "#880000", "#0000FF", "#006400", "#00FF00"), lwd=2),
                  title_position = "topleft", title = "Copy Number\n(State)")

  svs = Legend(labels = c('foldback', 'unbalanced', 'duplication',
                          'deletion', 'inversion', 'balanced'), type = "points",
               legend_gp = gpar(col =c(2, 3, 4, 1, 6, 8), lwd = 2), title_position = "topleft",
               title = "Structural Variants")

  draw(cn_leg, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  draw(svs, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"),
       just = c("right", "bottom"))


  #draw(gene_annotations,x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "top"))


}

#initialize a circos plot and fill with copy number ring, svs and legend
#copy_number -> copy number data.frame
#annotations -> gene annotations data.frame
#svs -> sv data.frame
circos_plot = function(copy_number, annotations, svs, label){
  require(circlize)

  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22, "X")))

  plot_copy_number(copy_number, annotations)
  plot_svs(svs)
  circos.clear()

  add_legend()
  title(label)
}

#make a pdf and write circos plot to it
#copy_number -> copy number data.frame
#annotations -> gene annotations data.frame
#svs -> sv data.frame
#outfile -> path for output pdf file
make_circos = function(copy_number, annotations, svs, outfile, sample_label){

  #save plot to pdf
  #pdf(outfile, 8, 9)

  #create plot
  circos_plot(copy_number, annotations, svs, sample_label)

  #turn off pdf
  #dev.off()

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
