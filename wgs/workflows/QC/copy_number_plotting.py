import pandas as pd
from matplotlib.lines import Line2D
import numpy as np


def read(copy_number):
    '''
    read in copy number data
    '''
    return pd.read_csv(copy_number, sep="\t")


def plot_anno_genes(anno_genes, min, max, ax):
    '''
    plot gene annotations
    '''
    for gene, gene_info in anno_genes.iterrows():

        x = [gene_info["pos"]/1000000, gene_info["pos"]/1000000]
        y = [min, max]
        ax.plot(x, y, color=gene_info["color"], zorder=1)

    return ax


def get_gene_annotation_data(chrom):
    '''
    make a hard coded dataframe of gene annotations
    '''
    anno_genes = {"CCNE1": [19, 30315215, "blue"],         "ERBB2": [17, 37886679, "green"],
                  "KRAS": [12, 25403870, "red"],           "MYC": [8, 128753674, "cyan"],
                  "PIK3CA": [3, 178957881, "magenta"],     "MECOM": [3, 169381406, "yellow"],
                  "RB1": [13, 49056122, "lightsteelblue"], "PTEN": [10, 89731687, "tan"],
                  "BRCA1": [17, 41277500, "black"],        "BRCA2": [13, 32973805, "grey"],
                  "RAD51C": [17, 56811703, "saddlebrown"], "PALB2": [16, 23652631, "pink"]}

    anno_genes = pd.DataFrame(anno_genes).T
    anno_genes = anno_genes.rename(columns={0: "chrom", 1: "pos", 2: "color"})
    return anno_genes[anno_genes.chrom == chrom]


def bin(positions, copy_number, state, n_bins, start, extent):
    '''
    bin coverage data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(positions, bins)

    position = [positions[digitized == i].mean() for i in range(1, len(bins))]
    lr = [copy_number[positions[digitized == i].index].mean() for i in range(1, len(bins))]
    state = [state[positions[digitized == i].index].mean() for i in range(1, len(bins))]
    return pd.DataFrame({"Position": position,
                       "LogRatio": lr,
                        "state": state})


def prepare_at_chrom(copy_number, chrom):
    '''
    prep copy number rdata for plotting at a chrom
    '''
    return copy_number[copy_number["Chr"] == chrom][["Position", "LogRatio"]]


def make_annotation_legend(anno_genes, axis):
    '''
    plot a legend for gene annotations
    '''
    colors = anno_genes.color.tolist()
    labels = anno_genes.index.tolist()
    lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors]
    l4 = axis.legend(lines, labels, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                    mode="expand", borderaxespad=0, ncol=6)
    bbox =axis.get_position()
    offset = -.04
    axis.set_position([bbox.x0, bbox.y0 + offset, bbox.x1-bbox.x0, bbox.y1 - bbox.y0])
    return axis


def plot(prepped_copy_number, anno_genes, axis):
    '''
    plot copy number data on axis
    '''
    axis.scatter(prepped_copy_number.Position / 1000000,
                 prepped_copy_number.LogRatio, s=0.1, color='black', zorder=2)

    axis.set_ylabel("CN")
    axis.set_ylim(-4, 6)

    if not anno_genes.empty:
        axis = plot_anno_genes(anno_genes, *axis.get_ylim(), axis)
        axis = make_annotation_legend(anno_genes, axis)

    return axis


def make_for_circos(copy_number, outfile):
    '''
    prep copy number data to be used by circos.r
    :param copy_number: input copy_number file
    :param outfile: path to prepped copy_number csv
    :return:
    '''
    df = read(copy_number)
    output = []
    chroms = df.Chr.unique()

    for chrom in chroms:

        prepped = copy_number[copy_number["Chr"] == chrom]
        prepped = bin(prepped.Position, prepped.LogRatio, prepped.TITANstate, 200,
                      prepped.Position.min(), prepped.Position.max())
        prepped["Chr"] = [chrom] * len(prepped.index)
        output.append(prepped)

    output = pd.concat(output)
    #the NaNs over centromere breaks the R code
    output = output[~output.Position.isna()]
    output.to_csv(outfile, index=False, header=True, sep="\t")


#EXAMPLE RUN
#run a titan_markers file through
# make_for_circos() to get input for circos.r