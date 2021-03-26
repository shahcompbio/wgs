import numpy as np
import pandas as pd
import pypeliner
import pysam
from wgs.workflows.sample_qc import genome_wide_plot
from wgs_qc_utils.plotter import gene_annotation_plotting
from wgs_qc_utils.reader import read_titan, read_remixt
import gzip
from wgs.utils import helpers
import os

def get_gene_annotations(outfile):
    chroms = list(map(str, range(1, 22))) + ["X"]
    annotations = pd.concat([gene_annotation_plotting.get_gene_annotation_data(chrom) for chrom in chroms])
    annotations.to_csv(outfile, sep="\t", index=False)


def circos(titan_calls, remixt_calls, sample_id, sv_calls,
           circos_plot_remixt, circos_plot_titan, tempdir, docker_image=None):

    helpers.makedirs(tempdir)

    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'scripts','circos.R')

    cmd = ['Rscript',script_path, titan_calls, remixt_calls, sv_calls,
           circos_plot_remixt, circos_plot_titan, sample_id]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def prep_data_for_circos(titan, remixt, sample_id, prepped_titan_calls,
                         prepped_remixt_calls):
    read_titan.make_for_circos(titan, prepped_titan_calls)

    read_remixt.make_for_circos(remixt, sample_id, prepped_remixt_calls)


def bin_data(positions, copy_number, state, n_bins, start, extent):
    '''
    bin coverage data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(positions, bins)

    position = [positions[digitized == i].mean() for i in range(1, len(bins))]
    lr = [copy_number[positions[digitized == i].index].mean() for i in range(1, len(bins))]
    state = [state[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame(
        {"Position": position,
         "LogRatio": lr,
         "state": state}
    )


def generate_coverage_bed(ref, bins_out, chromosomes, bins_per_chrom=2000):
    fasta = pysam.FastaFile(ref)

    chroms = dict(zip(fasta.references, fasta.lengths))
    chroms = {k: v for k, v in chroms.items() if k in chromosomes}

    chroms_all = []
    starts_all = []
    ends_all = []

    for chrom, length in chroms.items():
        step_size = int(length / bins_per_chrom)

        starts = [str(int(i * step_size) + 1) for i in range(bins_per_chrom)]
        ends = [str(int((i + 1) * step_size)) for i in range(bins_per_chrom)]

        assert len(starts) == len(ends)
        chroms = [chrom] * len(starts)

        starts_all += starts
        ends_all += ends
        chroms_all += chroms

    out = pd.DataFrame({"chrom": chroms_all, "starts": starts_all, "ends": ends_all})

    out = out[out.chrom == chromosomes]

    out.to_csv(bins_out, sep="\t", index=False, header=False)


def prep_sv_for_circos(sv_calls, outfile):
    svs = pd.read_csv(sv_calls, sep=",", dtype={'chromosome_1': str, 'chromosome_2': str})

    svs = svs[[
        'chromosome_1', 'position_1', 'chromosome_2', 'position_2',
        'rearrangement_type']
    ]

    types = [
        'foldback', 'unbalanced', 'duplication',
        'deletion', 'inversion', 'balanced'
    ]

    colors = [2, 3, 4, 1, 6, 8]
    svs["color"] = svs.rearrangement_type.replace(types, colors)

    svs.to_csv(outfile, index=False, header=True, sep="\t")


def roh_needs_parse(roh_calls):

    if roh_calls.endswith(".gz"):
        header = next(gzip.open(roh_calls, "rt"))
    else:
        header = next(open(roh_calls))

    if header.startswith("#"):
        return True
    
    return False

def parse_roh(roh_calls, parsed): 
    if roh_needs_parse(roh_calls):
        lines = [l for l in open(roh_calls) if "ST" in l and "#" not in l]
        with open(parsed, 'w') as f:
            f.write("%s\n" % "type,sample,chromosome,start,state,quality")
            for line in lines: 
                line = line.replace("\t", ",")
                f.write("%s" % line)
        f.close()



def samtools_coverage(
        bam_file, bed_file, output,
        mapping_qual, docker_image=None
):
    command = [
        "samtools", "bedcov", bed_file, bam_file,
        "-Q", mapping_qual, ">", output
    ]

    pypeliner.commandline.execute(*command, docker_image=docker_image)

    df_out = pd.read_csv(output, sep="\t", names=["chrom", "start", "end", "sum_cov"])
    df_out.to_csv(output, sep="\t", index=False, header=True)


def clear_header_label(f):
    data = pd.read_csv(f, sep="\t")
    data.to_csv(f, sep="\t", index=False, header=False)


def genome_wide(
        sample_id, titan, roh, germline_calls, somatic_calls, remixt,
        tumour_coverage, normal_coverage, breakpoints, chromosomes, pdf, normal_only=False
):
    if normal_only:
        genome_wide_plot.genome_wide_plot(
            None, sample_id, None, roh, germline_calls, None,
            None, normal_coverage, None, chromosomes, pdf, normal_only=normal_only
        )
    else:
        genome_wide_plot.genome_wide_plot(
            remixt, sample_id, titan, roh, germline_calls, somatic_calls,
            tumour_coverage, normal_coverage, breakpoints, chromosomes, pdf, normal_only=normal_only
        )
