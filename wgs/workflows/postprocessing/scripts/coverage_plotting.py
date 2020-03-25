import numpy as np
import pandas as pd


def prepare_at_chrom(coverage, chrom, bin=False):
    '''
    prepare coverage data at a chromosome to be plotted
    '''
    coverage = coverage[coverage.chrom == str(chrom)]

    if bin:
        coverage = bin_coverages(coverage.start, coverage.coverage,
                                 2000, coverage.start.min(), coverage.start.max())

    high = coverage.coverage.quantile(0.85)
    low = coverage.coverage.quantile(0.05)
    # preserve NaNs so as not to plot over centromere
    return coverage[(coverage.coverage.between(low, high))
                    | (coverage.coverage.isnull())]


def plot(prepped_coverage, ylim_min, ylim_max, axis, name):
    '''
    plot coverage data on an axis
    '''
    axis.plot(prepped_coverage.start / 1000000, prepped_coverage.coverage, color="black")
    axis.set_ylabel(name)

    axis.set_ylim(ylim_min, ylim_max)

    return axis


def read(coverage):
    '''
    read in coverage data
    '''
    cov = pd.read_csv(coverage, na_values="nan", )

    cov = cov.astype({"chrom": str})
    return cov


def bin_coverages(positions, coverages, n_bins, start, extent):
    '''
    bin coverage data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(positions, bins)

    start = [positions[digitized == i].min() for i in range(1, len(bins))]
    end = [positions[digitized == i].max() for i in range(1, len(bins))]
    coverage = [coverages[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame({"start": start,
                         "end": end,
                         "coverage": coverage})
