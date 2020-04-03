import os

import numpy as np
import pandas as pd
import pypeliner
import pysam
from scripts import plot
from wgs.utils import helpers

from wgs.utils import helpers

def circos(cn_calls, sv_calls, circos_plot, tempdir, annotations, docker_image=None):
    helpers.makedirs(tempdir)

    prepped_cn_calls = os.path.join(tempdir, 'prepped_cn_calls.csv')
    prep_cn_for_circos(cn_calls, prepped_cn_calls)

    prepped_sv_calls = os.path.join(tempdir, 'prepped_sv_calls.csv')
    prep_sv_for_circos(sv_calls, prepped_sv_calls)

    cmd = ['circos.R', annotations, prepped_cn_calls, prepped_sv_calls, circos_plot]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def bin_data(positions, copy_number, state, n_bins, start, extent):
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


def prep_cn_for_circos(copy_number, outfile):
    '''
    prep copy number data to be used by circos.r
    :param copy_number: input copy_number file
    :param outfile: path to prepped copy_number csv
    :return:
    '''
    copy_number = pd.read_csv(copy_number, sep="\t", dtype={'Chr': str})

    output = []
    chroms = copy_number.Chr.unique()

    for chrom in chroms:
        prepped = copy_number[copy_number["Chr"] == chrom]
        prepped = bin_data(prepped.Position, prepped.LogRatio, prepped.TITANstate, 200,
                           prepped.Position.min(), prepped.Position.max())
        prepped["Chr"] = [chrom] * len(prepped.index)
        output.append(prepped)

    output = pd.concat(output)
    # the NaNs over centromere breaks the R code
    output = output[~output.Position.isna()]
    output.to_csv(outfile, index=False, header=True, sep="\t")


def prep_sv_for_circos(sv_calls, outfile):
    svs = pd.read_csv(sv_calls, sep=",", dtype={'chromosome_1': str, 'chromosome_2': str})

    svs = svs[['chromosome_1', 'position_1',
               'chromosome_2', 'position_2', 'rearrangement_type']]

    types = ['foldback', 'unbalanced', 'duplication',
             'deletion', 'inversion', 'balanced']
    colors = [2, 3, 4, 1, 6, 8]
    svs["color"] = svs.rearrangement_type.replace(types, colors)

    svs.to_csv(outfile, index=False, header=True, sep="\t")


def count_depth(regions, input):
    """
    Count the depth of the read. For each genomic coordonate return the
    number of reads
    -----
    Parameters :
        chr : (str) name of the chromosome
    -----
    Returns :
        none
    """
    bamfile = pysam.AlignmentFile(input, 'rb')

    output_values = []

    for region in regions:
        chrom, start, stop = region.split('_')

        start = int(start)
        stop = int(stop)

        size = stop - start

        running_count = 0

        for pileupcolumn in bamfile.pileup(chrom, start, stop, stepper='nofilter'):
            running_count += pileupcolumn.nsegments

        output_values.append((chrom, start, stop, running_count / size))

    return output_values


def samtools_coverage(bam_file, output, chromosomes, reference, bins_per_chrom=2000):
    intervals = generate_intervals(reference, chromosomes, bins_per_chrom=bins_per_chrom)

    counts = count_depth(intervals, bam_file)

    with helpers.GetFileHandle(output, 'w') as outfile:
        outfile.write('chrom,start,stop,coverage\n')
        for value in counts:
            value = ','.join(map(str, value)) + '\n'
            outfile.write(value)


def generate_intervals(ref, chromosomes, bins_per_chrom=2000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue

        step_size = int(length / bins_per_chrom)

        for i in range(bins_per_chrom):
            start = str(int(i * step_size) + 1)
            end = str(int((i + 1) * step_size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def genome_wide_plot(
        cn_markers, roh_calls, germline_vcf, somatic_vcf,
        tumour_coverage, normal_coverage, pdf_output
):
    plot.main(
        cn_markers, roh_calls, germline_vcf, somatic_vcf,
        tumour_coverage, normal_coverage, pdf_output
    )
