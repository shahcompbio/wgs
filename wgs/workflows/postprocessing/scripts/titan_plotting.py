import pandas as pd
import numpy as np
import gene_annotation_plotting
from matplotlib.lines import Line2D

def read(copy_number):
    """
    read in copy number data
    :param copy_number: copy_number input data (titan_markers.csv.gz)
    :returns: read-in copynumber pandas dataframe with plot colors
    """

    read = pd.read_csv(copy_number, sep="\t")
    read = read.astype({"Chr": str})

    n_extra = read.TITANstate.max() - 7

    colors = ["#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000"] + ["#FF0000"] * n_extra

    read["color"] = read.TITANstate.apply(lambda state: colors[state])

    return read


def prepare_at_chrom(copy_number, chrom):
    """
    get copy number data at a chromosome
    :param copy_number: read in copy number data
    :param chrom: chromosome to extract (str)
    :return: copy number data at chrom
    """
    '''
    prep copy number rdata for plotting at a chrom
    '''

    if isinstance(chrom, int):
        chrom = str(chrom)

    return copy_number[copy_number["Chr"] == chrom][["Position", "LogRatio", "color"]]

def add_titan_legend(axis):
    labels = ["1", "2", "3", "4", "5",
              "6", "7", "8", "9+"]

    colors = ["#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", "#FF0000"]

    lines = [Line2D([0], [0], marker='o', color='w', label='major axis',
                    markerfacecolor=c, markersize=10) for c, _ in zip(colors, labels)]

    axis.legend(lines, labels, ncol=3, columnspacing=0.02,
                loc="center left", title="Titan", frameon=False,
                borderpad=0, borderaxespad=0)

    return axis


def plot(prepped_copy_number, anno_genes, axis, chrom_max):
    """
    plot prepped copy number data on axis
    :param prepped_copy_number: prepped copy number data (read->prepare_at_chrom->
    :param anno_genes:
    :param axis:
    :return:
    """
    '''
    plot copy number data on axis
    '''
    axis.set_ylim(-4, 6)
    axis.set_xlim(0, chrom_max)
    axis.grid(True, linestyle=':')
    axis.spines['left'].set_position(('outward', 5))
    axis.spines['bottom'].set_position(('outward', 5))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_yticks(np.arange(-4, 6, 1))
    axis.set_xticks(np.arange(0, prepped_copy_number.Position.max() / 1000000, 25))

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.scatter(prepped_copy_number.Position / 1000000,
                 prepped_copy_number.LogRatio, s=0.1, color=prepped_copy_number.color)

    axis.set_ylabel("Titan", fontsize=14, fontname="Arial")

    if not anno_genes.empty:
        axis = gene_annotation_plotting.plot_anno_genes(anno_genes, *axis.get_ylim(), axis)


    return axis


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for circos plotting ~~~~~~~~~~~~~~~~~~~ #


def bin_copy_number(copy_number, n_bins=200):
    """
    bin copy number data
    :param positions:
    :param copy_number:
    :param state:
    :param n_bins:
    :param start:
    :param extent:
    :return:
    """
    positions = copy_number.Position
    state = copy_number.TITANstate
    log_ratio = copy_number.LogRatio

    bins = np.linspace(copy_number.Position.min(), copy_number.Position.max(), n_bins)
    digitized = np.digitize(positions, bins)

    position = [positions[digitized == i].mean() for i in range(1, len(bins))]
    lr = [log_ratio[positions[digitized == i].index].mean() for i in range(1, len(bins))]
    state = [state[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame({"Position": position, "LogRatio": lr,
                         "state": state})


def make_for_circos(copy_number, outfile):
    """
    prepare cn data for circos plot
    bin the data at each chromosome, get rid of nans + unused cols
    :param copy_number: copy number input (titan_markers.csv.gz)
    :param outfile: output path for prepped data
    :return: prepped cn data ready to go for circos plot
    """

    cn = read(copy_number)

    chroms = cn.Chr.unique()
    output = dict(zip(chroms, [[]] * len(chroms)))

    for chrom in chroms:
        cn_at_chrom = cn[cn.Chr == chrom]
        out = bin_copy_number(cn_at_chrom, n_bins=200)

        out["Chr"] = [chrom] * len(out.index)
        output[chrom] = out

    output = pd.concat(output)
    output = output[~output.Position.isna()]

    output.to_csv(outfile, index=False, header=True, sep="\t")