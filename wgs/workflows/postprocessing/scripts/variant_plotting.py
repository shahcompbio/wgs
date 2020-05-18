import pandas as pd
import numpy as np
import gzip


def plot_bar(variants, axis, name, chrom_max):
    '''
    plot variants on axis
    '''
    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.bar(variants.location / 1000000,  variants.n_events,
             facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    return axis


def plot_fill(variants, axis, name, chrom_max):
    '''
    plot variants on axis
    '''
    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False


    axis.fill_between(variants.location / 1000000, 0, variants.n_events,
                      facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    return axis


def read_consensus(f):
    f = pd.read_csv(f)
    f = f[["chrom", "pos", "AC_NORMAL", "DP_NORMAL", "AC_TUMOUR", "DP_TUMOUR"]]
    f["VAF_normal"] = f.AC_NORMAL/f.DP_NORMAL
    f["VAF_tumor"] = f.AC_TUMOUR/f.DP_TUMOUR
    f = f.drop(["AC_NORMAL", "DP_NORMAL", "AC_TUMOUR", "DP_TUMOUR"], axis=1)
    f = f.rename(columns = {"chrom":"chr"})
    f = f.astype({"chr":str})
    return f


def read_titan_vcf(f):
    data = pd.DataFrame(parse(open(f), "\t", 11), columns=["chr", "pos", "id", "ref", "alt",
                                                          "qual", "filter", "info",
                                                           "format", "tumour", "normal"])

    data = data.astype({"pos": np.int64, "chr": str})[["chr", "pos","format",
                                                           "tumour", "normal"]]
    cols = data.format.str.split(":").tolist()[0]

    cols_norm = [c + "_norm" for c in cols]
    data[cols_norm] = data.normal.str.split(":", expand=True)
    cols_tum= [c + "_tum" for c in cols]
    data[cols_tum] = data.tumour.str.split(":", expand=True)

    data = data.drop("format", axis=1)
    data = data.drop("normal", axis=1)

    data[["allele_1", "allele_2"]] = data.GT_norm.str.split("/", expand=True)

    return data


def prepare_at_chrom(variants, chrom, n_bins=200):
    '''
    prepare variants data to be plotted at a chrom
    '''

    variants = variants[variants["chr"] == str(chrom)]
    return bin_frequencies(variants.pos, n_bins, variants.pos.min(),
                           variants.pos.max())


def bin_frequencies(locations, n_bins, start, extent):
    '''
    bin variant data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(locations, bins)

    binned_loc = [locations[digitized == i].mean() for i in range(1, len(bins))]
    n_events = [len(locations[digitized == i]) for i in range(1, len(bins))]

    return pd.DataFrame({"location": binned_loc,
                       "n_events": n_events})

def read_full(f):

    data = pd.DataFrame(parse(open(f), "\t", 10), columns=["chr", "pos", "id", "ref", "alt",
                                                          "qual", "filter", "info", "format",
                                                             "normal"])
    data = data.astype({"pos": np.int64, "chr": str})[["chr", "pos","format",
                                                           "normal"]]
    cols = data.format.str.split(":").tolist()[0]

    cols_norm = [c + "_norm" for c in cols]
    data[cols_norm] = data.normal.str.split(":", expand=True)
    for c in cols_norm:
        if c != "GT_norm":
            data = data.drop(c, axis=1)

    data = data.drop("format", axis=1)
    data = data.drop("normal", axis=1)

    data[["allele_1", "allele_2"]] = data.GT_norm.str.split("/", expand=True)
    return data


def read(f):
    '''
    read in
    '''

    data = pd.DataFrame(parse(f, "\t", 2), columns=["chr", "pos"])
    data = data.astype({"pos": np.int64, "chr": str})
    return data[["chr", "pos"]]


def parse(f, sep, n):
    """
    parse vcf
    :param f: vcf file
    :param sep: seperator to use
    :param n: number of elements to take per line
    :return: parsed vcf file as list of lists
    """
    if f.endswith(".gz"):
        with gzip.open(f, 'rt') as f:
            return [line.split(sep)[:n] for line in f
                    if not line.startswith('#')]
    else:
        return [line.split(sep)[:n] for line in open(f)
                if not line.startswith('#')]
