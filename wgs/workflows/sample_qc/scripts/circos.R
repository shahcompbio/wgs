library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

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

#plot remixt
plot_remixt = function(copy_number){
  print('plot_remixt 1')
  copy_number$chrom = gsub("chr", "", copy_number$chrom) # chr2 -> 2
  chrom_maxes = rep(-1, length(unique(copy_number$chrom)))
  names(chrom_maxes) = unique(copy_number$chrom)
  for (chrom in names(chrom_maxes)){
    print('plot_remixt: chrom loop')
    print(chrom)
    chrom_maxes[chrom] = max(copy_number[copy_number$chrom==chrom, ]$start)
    #chrom_maxes[chrom] = max(copy_number[copy_number$chrom==chrom, ]$start)
  }
  print('plot_remixt 2')
  #chrom_starts = rep(0, length(chroms_needed)) #
  chrom_starts = rep(0, length(chrom_maxes)) #
  print('names(chrom_starts)')
  print(names(chrom_starts))
  print('chrom_maxes')
  print(chrom_maxes)
  print('names(chrom_maxes)')
  print(names(chrom_maxes))
  print('length(chrom_maxes)')
  print(length(chrom_maxes))
  print('chrom_maxes')
  print(chrom_maxes)
  names(chrom_starts) = names(chrom_maxes)
  print('names(chrom_starts) 2')
  print(names(chrom_starts))
  print('length(chrom_starts) 2')
  print(length(chrom_starts))
  print('chrom_starts 2')
  print(chrom_starts)
  chrom_positions = c(chrom_starts, chrom_maxes)
  print('plot_remixt 3')

  circos.track(factors=paste0("chr", copy_number$chrom),
               ylim=c(0, 8), "track.height" = 0.15)

  print('plot_remixt 4')
  print('chrom_positions')
  print(chrom_positions)
  chrom_factors = paste0("chr",names(chrom_positions))
  print('chrom_factors')
  print(chrom_factors)
  print('length(chrom_factors)')
  print(length(chrom_factors))
  print('unname(chrom_positions)')
  print(unname(chrom_positions))
  print('length(unname(chrom_positions))')
  print(length(unname(chrom_positions)))
  for (state in 0:8){
    circos.trackLines(factors = chrom_factors,
                      x = unname(chrom_positions), 
                      #y = rep(state, 46), 
                      y = rep(state, length(chrom_positions)), 
                      col="grey")

  }
  print('plot_remixt 5')


  cn_combined = copy_number
  cn_combined["combined"] = cn_combined$major_raw_e + cn_combined$minor_raw_e
  cn_combined = cn_combined[cn_combined$combined <= 8,]

  circos.trackPoints(factors= paste0("chr", cn_combined$chrom),
                     x = cn_combined$start, y = cn_combined$combined, cex = 0.1, col = "black")

  circos.track(factors= paste0("chr", cn_combined$chrom),
               ylim=c(0, 8), "track.height" = 0.15)
  print('plot_remixt 6')
  for (state in 0:8){
    circos.trackLines(factors = paste0("chr",names(chrom_positions)),
                      x = unname(chrom_positions), 
                      #y = rep(state, 46), 
                      y = rep(state, length(chrom_positions)), #
                      col="grey")

  }


  print('plot_remixt 7')
  circos.trackPoints(factors= paste0("chr", copy_number$chrom),
                     x = copy_number$start, y = copy_number$major_raw_e, cex = 0.1, col = "red")

  circos.trackPoints(factors= paste0("chr", copy_number$chrom),
                     x = copy_number$start, y = copy_number$minor_raw_e, cex = 0.1, col = "blue")
}


plot_titan = function(copy_number){
  print('plot_titan 1')
  cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
  names(cnCol) <- c(0:500)

  copy_number$Chr = gsub("chr", "", copy_number$Chr)
  chroms_needed = c(1:22, "X", "Y")
  if (! all (chroms_needed %in% copy_number$Chr)){
    copy_number = fill_cn_data_chroms(copy_number, "Chr", chroms_needed)
    print('plot_titan: copy_number')
    print(copy_number)
  }

  chrom_maxes = rep(-1, length(unique(copy_number$Chr)))
  names(chrom_maxes) = unique(copy_number$Chr)
  print('plot_titan: chrom_maxes')
  print(chrom_maxes)
  for (chrom in names(chrom_maxes)){
    chrom_maxes[chrom] = max(copy_number[copy_number$Chr==chrom, ]$Position)
  }
  print('plot_titan: chrom_maxes 2')
  print(chrom_maxes)
  chrom_starts = rep(0, length(chroms_needed))

  names(chrom_starts) = names(chrom_maxes)
  print('plot_titan: chrom_starts')
  print(chrom_starts)
  chrom_positions = c(chrom_starts, chrom_maxes)
  print('plot_titan: chrom_positions')
  print(chrom_positions)

  circos.track(factors= paste0("chr", copy_number$Chr),
               ylim=c(0, 8), "track.height" = 0.15)
  print('plot_titan 2')

  for (state in 0:8){
    circos.trackLines(factors = paste0("chr",names(chrom_positions)),
                      x = unname(chrom_positions), 
                      y = rep(state, length(chrom_positions)), 
                      col="grey")
  }
  print('plot_titan 3')

  for (state in unique(copy_number$state)){
      x = copy_number$Position[copy_number$state == state]
      y = copy_number$LogRatio[copy_number$state == state]
      f = copy_number$Chr[copy_number$state == state]
      col = cnCol[state+1]
      circos.trackPoints(factors= paste0("chr", f),  x = x, y = y, cex = 0.1, col=col)
  }
  print('plot_titan 4')

}

#add the svs to the circos plot as links
#svs -> data.frame of breakpoints
plot_svs = function(svs, palette){
  if (! nrow(svs)){
    return()
  }
  svs$chromosome_1 = gsub('chr', '', svs$chromosome_1)
  svs$chromosome_2 = gsub('chr', '', svs$chromosome_2)
  for (sv in 1:nrow(svs)){
    chr1 = paste0("chr", svs[sv, "chromosome_1"])
    print('plot_svs: chr1')
    print(chr1)
    chr2 = paste0("chr", svs[sv, "chromosome_2"])
    print('plot_svs: chr2')
    print(chr2)
    p1 =  svs[sv, "position_1"]
    print('plot_svs: p1')
    print(p1)
    p2 =  svs[sv, "position_2"]
    print('plot_svs: p2')
    print(p2)
    color = palette[svs[sv, "rearrangement_type"]]
    print('plot_svs: color')
    print(color)

    if (chr1 == "chrY" | chr2 == "chrY"){
      next
    }

    circos.link(chr1, p1, chr2, p2, col = color)

  }
}

#add a legend to the plot
add_legend = function(palette, remixt, titan){

  if (titan){
    titan_colors =  c("#00FF00", "#006400", "#0000FF", "#880000",
                      "#BB0000", "#CC0000", "#DD0000", "#EE0000", "#FF0000")
    names(titan_colors) = c("1", "2", "3", "4", "5", "6", "7", "8", "9+")

    cn_leg = Legend(labels = rev(names(titan_colors)), type = "points",
                 legend_gp = gpar(col = rev(titan_colors), lwd = 2), title_position = "topleft",
                 title = "Titan State")


    cn_states = Legend(labels = c("Titan LogRatio (-4 - 6)"), type="lines",
                       legend_gp=gpar(col=c("gray"), lwd=2),
                       title_position = "topleft")
  }

  if (remixt){
    cn_leg = Legend(labels = c("major + minor","major", "minor"), type="points",
                  legend_gp=gpar(col=c("black","red", "blue"), lwd=2),
                  title_position = "topleft", title = "Remixt")

    cn_states = Legend(labels = c("Copy Number State (0 - 8+)"), type="lines",
                       legend_gp=gpar(col=c("gray"), lwd=2),
                       title_position = "topleft")
  }


  cn_leg = packLegend(cn_leg, cn_states)


  svs = Legend(labels = names(palette), type = "points",
               legend_gp = gpar(col = palette, lwd = 2), title_position = "topleft",
               title = "Structural Variants")

  draw(cn_leg, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  draw(svs, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"),
       just = c("right", "bottom"))

}

#initialize a circos plot and fill with copy number ring, svs and legend
#copy_number -> copy number data.frame
#annotations -> gene annotations data.frame
#svs -> sv data.frame
circos_plot_titan = function(titan, svs, label){
  require(circlize)
  print('circos_plot_titan 1')

  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  print('circos_plot_titan 2')
  palette = brewer.pal(6, "Set2")
  names(palette) = c("duplication", "deletion", "unbalanced", "balanced", "foldback", "inversion")

  plot_titan(titan)
  plot_svs(svs, palette)
  circos.clear()
  par(family="sans")
  add_legend(palette, FALSE, TRUE)
  title(label)
  print('circos_plot_titan 3')
}


circos_plot_remixt = function(remixt, svs, label){
  print('circos_plot_remixt 1')
  require(circlize)

  print('circos_plot_remixt 2')
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  print('circos_plot_remixt 3')
  palette = brewer.pal(6, "Set2")
  names(palette) = c("duplication", "deletion", "unbalanced", "balanced", "foldback", "inversion")

  print('circos_plot_remixt 4')
  plot_remixt(remixt)
  print('circos_plot_remixt 5')
  plot_svs(svs, palette)
  print('circos_plot_remixt 6')
  circos.clear()
  print('circos_plot_remixt 7')
  par(family="sans")
  print('circos_plot_remixt 8')
  add_legend(palette, TRUE, FALSE)
  print('circos_plot_remixt 9')
  title(label)
  print('circos_plot_remixt 10')
}



make_circos_remixt = function(copy_number, svs, outfile, sample_label){
  print('make_circos_remixt 1')

  #save plot to pdf
  pdf(outfile, 8, 9)
  print('make_circos_remixt 2')

  #create plot
  circos_plot_remixt(copy_number, svs, sample_label)
  print('make_circos_remixt 3')

  #turn off pdf
  dev.off()
  print('make_circos_remixt 4')

}

make_circos_titan = function(copy_number, svs, outfile, sample_label){

  #save plot to pdf
  pdf(outfile, 8, 9)

  #create plot
  circos_plot_titan(copy_number, svs, sample_label)

  #turn off pdf
  dev.off()

}


args <- commandArgs(TRUE)
titan <- read.csv(args[1], sep = "\t")
titan$Chr <- gsub('x', 'X', titan$Chr)
titan$Chr <- gsub('y', 'Y', titan$Chr)

remixt <- args[2]

print('main:remixt') # prints out a path str
print(remixt) #
if (remixt != "NULL"){
    remixt <- read.csv(args[2], sep = "\t")
    remixt$chrom <- gsub('x', 'X', remixt$chrom)
    remixt$chrom <- gsub('y', 'Y', remixt$chrom)
}

sv_calls <- read.csv(args[3], colClasses=c("rearrangement_type"="character"))
remixt_name <- args[4]
titan_name <- args[5]

sample_id <- args[6]


print('main:remixt 2') # prints out a table - no problem until here
if (remixt != "NULL"){
    print('main:just before make_circos_remixt')
    make_circos_remixt(remixt, sv_calls, remixt_name, sample_id)

}
make_circos_titan(titan, sv_calls, titan_name, sample_id)
