import pandas as pd
from matplotlib.lines import Line2D
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ code taken from Andrew McPherson ~~~~~~~~~~~~~~~~~~~~~~ #

def read(remixt_file, sample_label):
    """
    read in remixt data into pandas dataframe
    :param remixt_file: remixt input file (i.e. results_files_*.h5)
    :param sample_label: label of sample inside h5
    :return: parsed pandas dataframe
    """
    cnv_data = list()
    with pd.HDFStore(remixt_file, 'r') as store:
        prefix = '/copy_number/{pred_tool}/sample_{sample_id}'.format(sample_id=sample_label,
                                                                      pred_tool=sample_label)
        mixture = store['/mix']
        cn = store['/cn']
        # Ensure that major and minor are correct, otherwise swap
        swap = (cn['major_1'] < cn['minor_1']) & (cn['major_2'] < cn['minor_2'])
        cn.loc[swap, ['major_1', 'minor_1']] = cn.loc[swap, ['minor_1', 'major_1']].values
        cn.loc[swap, ['major_2', 'minor_2']] = cn.loc[swap, ['minor_2', 'major_2']].values
        cn.loc[swap, 'major_is_allele_a'] = 1 - cn.loc[swap, 'major_is_allele_a']
        cn['total_1'] = cn['minor_1'] + cn['major_1']
        cn['total_2'] = cn['minor_2'] + cn['major_2']
        if mixture[1] > mixture[2]:
            cn['major'] = cn['major_1']
            cn['minor'] = cn['minor_1']
            cn['total'] = cn['total_1']
        else:
            cn['major'] = cn['major_2']
            cn['minor'] = cn['minor_2']
            cn['total'] = cn['total_2']
        swap = (cn['major'] < cn['minor'])
        cn.loc[swap, ['major', 'minor']] = cn.loc[swap, ['minor', 'major']].values
        assert (cn['major'] >= cn['minor']).all()
        cn['is_subclonal'] = (
                                     (cn['minor_1'] != cn['minor_2']) |
                                     (cn['major_1'] != cn['major_2'])) * 1
        cn['tumour_content'] = mixture[1:].sum()
        cnv_data.append(cn)
    cnv_data = pd.concat(cnv_data, ignore_index=True)
    cnv_data["total_raw_e"] = cnv_data.major_raw_e + cnv_data.minor_raw_e
    return cnv_data


def prepare_at_chrom(parsed_remixt, chrom):
    """
    prep copy number rdata for plotting at a chrom
    :param remixt: parsed remixt pandas dataframe
    :param chrom: chromosome
    :return: remixt at chromosome
    """
    return parsed_remixt[parsed_remixt["chromosome"] == chrom]


def add_remixt_legend(axis):
    labels = ["major axis", "minor axis", "snv copy number"]

    lines = [Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='r', markersize=10, alpha=0.4),
             Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='b', markersize=10, alpha=0.4),
             Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='black', markersize=5, alpha=1)]

    axis.legend(lines, labels, ncol=1, loc="center left", title="Remixt",
                frameon=False, borderpad=0, borderaxespad=0)

    return axis


def plot(prepped_remixt, axis, chrom_max):
    """
    plot remixt data on axis
    :param prepped_remixt: prepped remixt data
    :param anno_genes: gene annotations to add
    :param axis: axis to plot on
    :return: axis with plot added
    """
    axis.set_ylim(0, 8)
    axis.set_xlim(0, chrom_max)

    axis.grid(True, linestyle=':')
    axis.spines['left'].set_position(('outward', 5))
    axis.spines['bottom'].set_position(('outward', 5))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_yticks(range(0, 9))
    axis.set_xticks(np.arange(0, prepped_remixt.start.max() / 1000000, 25))

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    # return axis

    axis.scatter(prepped_remixt.start / 1000000,
                 prepped_remixt.major_raw,  c="red", s=20, marker="o", alpha=0.4)
    axis.scatter(prepped_remixt.start / 1000000,
                 prepped_remixt.minor_raw,  c="blue", s=20, marker="o", alpha=0.4)

    axis.set_ylabel("Remixt", fontsize=14, fontname="Arial")

    return axis


def make_for_circos(remixt, sample_id, prepped_remixt):
    remixt = read(remixt, sample_id)
    remixt.to_csv(prepped_remixt, sep="\t", index=False, header=True)

# import matplotlib.pyplot as plt
# f, ax = plt.subplots()
#
#
# a = read("/Users/abramsd/work/CODE/wgs/wgs/workflows/postprocessing/results_files_007.h5", "Sample_007")
# a = prepare_at_chrom(a, "1")
# plot(a,"a", ax)
# add_remixt_legend(ax)
# plt.show()