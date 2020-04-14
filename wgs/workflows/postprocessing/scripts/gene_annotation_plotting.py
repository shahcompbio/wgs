import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.legend import Legend


def add_gene_annotation_legend(anno_genes, axis):
    """
    add legend for gene annotations
    :param anno_genes: list of gene annotations to add
    :param axis: axis to plot on
    :return: axis with legend
    """

    colors = anno_genes.color.tolist()
    labels = anno_genes.index.tolist()
    lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors]


    axis.legend( lines, labels,
                 ncol=1, frameon=False, loc="center left", borderpad=0, borderaxespad=0, title="Gene Annotations")

    # axis.add_artist(leg)

    return axis


def plot_anno_genes(anno_genes, min, max, ax):
    """
    plot annotation genes on axis as vertical lines
    :param anno_genes: anno genes to plot
    :param min: minimum of y -axis
    :param max: maximum of y axis
    :param ax: axis to plot on
    :return: axis with added annotations
    """

    for gene, gene_info in anno_genes.iterrows():

        x = [gene_info["pos"]/1000000, gene_info["pos"]/1000000]
        y = [min, max]
        ax.plot(x, y, color=gene_info["color"])

    return ax


def get_gene_annotation_data(chrom):
    """
    gene gene annotation data at a chrom
    :param chrom: chromosome
    :return: gene annotations as pandas dataframe
    """

    anno_genes = {"CCNE1": [19, 30315215, "blue"], "ERBB2": [17, 37886679, "green"],
                  "KRAS": [12, 25403870, "red"], "MYC": [8, 128753674, "cyan"],
                  "PIK3CA": [3, 178957881, "magenta"], "MECOM": [3, 169381406, "yellow"],
                  "RB1": [13, 49056122, "lightsteelblue"], "PTEN": [10, 89731687, "tan"],
                  "BRCA1": [17, 41277500, "black"], "BRCA2": [13, 32973805, "grey"],
                  "RAD51C": [17, 56811703, "saddlebrown"], "PALB2": [16, 23652631, "pink"]}

    anno_genes = pd.DataFrame(anno_genes).T
    anno_genes = anno_genes.rename(columns={0: "chrom", 1: "pos", 2: "color"})
    anno_genes = anno_genes.astype({"chrom": str})

    return anno_genes[anno_genes.chrom == chrom]
