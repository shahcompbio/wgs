import os

import numpy as np
import pandas as pd
import pypeliner
import pysam

from .scripts import genome_wide_plot
from .scripts import titan_plotting
from .scripts import remixt_plotting
from .scripts import gene_annotation_plotting
from .scripts import coverage_plotting
from wgs.utils import helpers


def get_gene_annotations( outfile):

    chroms = list(map(str, range(1, 22))) + ["X"]
    annotations = pd.concat([gene_annotation_plotting.get_gene_annotation_data(chrom) for chrom in chroms])
    annotations.to_csv(outfile, sep="\t", index=False)


def circos(titan_calls, remixt_calls, sample_id, sv_calls,
           circos_plot_remixt, circos_plot_titan, tempdir, docker_image=None):

    helpers.makedirs(tempdir)

    ann_file = os.path.join(tempdir, "gene_annotations.tsv")
    get_gene_annotations(ann_file)

    prepped_titan_calls = os.path.join(tempdir, 'prepped_titan_calls.csv')
    titan_plotting.make_for_circos(titan_calls, prepped_titan_calls)

    prepped_remixt_calls = os.path.join(tempdir, 'prepped_remixt_calls.csv')
    remixt_plotting.make_for_circos(remixt_calls, sample_id, prepped_remixt_calls)

    cmd = ['circos.R', prepped_titan_calls, prepped_remixt_calls, ann_file, sv_calls,
           circos_plot_remixt, circos_plot_titan, sample_id]
    # cmd = ['ls']

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

    output_values = [None] * len(regions)

    for i, region in enumerate(regions):
        chrom, start, stop = region.split('_')

        start = int(start)
        stop = int(stop)
        size = stop - start

        running_count = 0

        for pileupcolumn in bamfile.pileup(chrom, start, stop, stepper='nofilter'):

            running_count += pileupcolumn.nsegments

        output_values[i] = (chrom, start, stop, running_count / size)

    return output_values


def samtools_coverage(bam_file, output, chromosomes, reference, bins_per_chrom=2000):
    # out = coverage_plotting.read("/juno/work/shah/abramsd/CODE/inputs/merged007_TT")
    # cov = coverage_plotting.prepare_at_chrom(out, chromosomes)
    #
    # cov.to_csv(output, sep="\t", index=False)
    intervals = generate_intervals(reference, chromosomes, bins_per_chrom=bins_per_chrom)

    counts = count_depth(intervals, bam_file)

    with helpers.GetFileHandle(output, 'w') as outfile:
        outfile.write('chrom,start,stop,coverage\n')
        for value in counts:
            value = ','.join(map(str, value)) + '\n'
            outfile.write(value)


def generate_intervals(ref, chromosome, bins_per_chrom=2000):
    fasta = pysam.FastaFile(ref)
    chrom_lengths = dict(zip(fasta.references, fasta.lengths))
    assert chromosome in chrom_lengths.keys()

    length = chrom_lengths[chromosome]

    intervals = [None] * bins_per_chrom

    step_size = int(length / bins_per_chrom)

    for i in range(bins_per_chrom):
        start = str(int(i * step_size) + 1)
        end = str(int((i + 1) * step_size))
        intervals[i] = str(chromosome) + "_" + start + "_" + end
    return intervals


def genome_wide(
        remixt, remixt_label, titan, roh, germline_calls, somatic_calls,
        tumour_coverage, normal_coverage, breakpoints, ideogram, chromosomes, pdf
):

    genome_wide_plot.genome_wide_plot(
        remixt, remixt_label, titan, roh, germline_calls, somatic_calls,
        tumour_coverage, normal_coverage, breakpoints, ideogram, chromosomes, pdf
    )
