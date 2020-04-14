import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def prepare_at_chrom(coverage, chrom, bin=False, n_bins=200):
    """
    prepare coverage data at a chromosome to be plotted
    :param coverage: read-in coverage data
    :param chrom: chromosome
    :param n_bins: number of bins to bin data into
    :return: binned dataframe
    """
    coverage = coverage[coverage.chrom == str(chrom)]

    if bin:
        coverage = bin_coverages(coverage, n_bins)

    high = coverage.coverage.quantile(0.99)
    low = coverage.coverage.quantile(0.01)

    # preserve NaNs so as not to plot over centromere
    return coverage[(coverage.coverage.between(low, high))
                    | (coverage.coverage.isnull())]


def plot(prepped_coverage, ylim_min, ylim_max, axis, name, chrom_max):
    """
    plot coverage data on an axis
    :param prepped_coverage: prepped coverage data
    :param ylim_min: min for y axis
    :param ylim_max: max for y axis
    :param axis: axis to plot on
    :param name: name for axis (i.e. normal or tumour coverage)
    :return: axis with plot
    # """
    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.fill_between(prepped_coverage.start / 1000000, ylim_min,
                      prepped_coverage.coverage, facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    axis.set_ylim(ylim_min, ylim_max)

    return axis


def read(coverage):
    """
    read in coverage data
    :param coverage: coverage data
    :return: pandas dataframe
    """

    cov = pd.read_csv(coverage, na_values="nan",
                           sep="\t")

    cov = cov.astype({"chrom": str})
    return cov


def bin_coverages(coverage, n_bins):
    """
    bin coverage data
    :param coverage: input coverage data
    :param n_bins: number of bins to seperate data into
    :return: binned coverage data as dataframe
    """
    positions = coverage.start
    coverages = coverage.coverage

    bins = np.linspace(positions.min(), positions.max(), n_bins)
    digitized = np.digitize(positions, bins)

    start = [positions[digitized == i].min() for i in range(1, len(bins))]
    end = [positions[digitized == i].max() for i in range(1, len(bins))]
    coverage = [coverages[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame({"start": start, "end": end, "coverage": coverage})

