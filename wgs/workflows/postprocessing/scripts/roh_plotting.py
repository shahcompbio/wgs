import numpy
import pandas as pd
from wgs.utils import helpers


def read(file):
    '''
    read in roh data
    '''

    skiprows = []

    with helpers.GetFileHandle(file) as infile:
        for i, line in enumerate(infile):
            if line.startswith('#') or line.startswith('RG'):
                skiprows.append(i)

    df = pd.read_csv(file, sep="\t", skiprows=skiprows)

    df.columns = ['ST', 'sample', 'chrom', 'pos', 'state', 'qual']
    df.chrom = df.chrom.astype('str')

    return df


def prepare_at_chrom(roh, chrom):
    '''
    prep roh data for plotting
    '''
    roh = roh[(roh.qual > 30) & (roh.chrom == chrom)]
    return bin(roh, 200, roh.pos.min(), roh.pos.max())


def bin(data, bins, start, end):
    '''
    bin roh data
    '''
    pos = data.pos
    state = data.state

    bins = numpy.linspace(start, end, bins)
    assigned = numpy.digitize(pos, bins)
    bin_range = range(1, len(bins))

    bin_start = [pos[assigned == i].min() for i in bin_range]

    bin_end = [pos[assigned == i].max() for i in bin_range]
    bin_mean = [pos[assigned == i].mean() for i in bin_range]

    # calculate the number of assigned hom states over total number and round it
    bin_state = [round(len(state[assigned == i][state[assigned == i] == 1])
                       / (len(state[assigned == i]) + 1)) for i in bin_range]

    n_homo_locs = [len(state[assigned == i][state[assigned == i] == 1])
                   for i in bin_range]

    n_het_locs = [len(state[assigned == i][state[assigned == i] == 0])
                  for i in bin_range]

    return pd.DataFrame({"bin_pos_mean": bin_mean, "bin_pos_start": bin_start,
                         "bin_pos_end": bin_end, "bin_overall_state": bin_state,
                         "n_homo_locs": n_homo_locs, "n_het_locs": n_het_locs})


def plot(roh, axis):
    '''
    plot roh data
    '''
    axis.plot(roh.bin_pos_mean / 1000000, roh.bin_overall_state, color='black')
    axis.set_ylabel("ROH")

    return axis

    # write_out_bin_roh_data(sample_36_roh, "bins/036.tsv")
