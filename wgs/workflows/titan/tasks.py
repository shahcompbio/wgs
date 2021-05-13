import os
import shutil
import tarfile

import numpy as np
import pandas as pd
import pypeliner
import pysam
from wgs.utils import csvutils
from wgs.utils import helpers
from wgs.utils import pdfutils
from wgs.utils import vcfutils
from wgs.workflows.titan.scripts import PygeneAnnotation
from wgs.workflows.titan.scripts import ReadCounter
from wgs.workflows.titan.scripts import parse_titan
from wgs.workflows.titan.scripts import vcf_to_counts


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size) + 1)
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def merge_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)


def convert_museq_vcf2counts(infile, outfile, het_positions):
    """
    Transform museq vcf file to a text file of counts

    :param infile: merged temporary vcf file from museq run (museq.vcf)
    :param outfile: temporary text file of counts (museq_postprocess.txt)
    """
    vcf_to_counts(infile, outfile, het_positions)


def run_readcounter(input_bam, output_wig, chromosomes, config):
    rc = ReadCounter(
        input_bam, output_wig, config['w'],
        chromosomes, config['q'],
    )
    rc.main()


def calc_correctreads_wig(tumour_wig, normal_wig, target_list, outfile, gc_wig, map_wig, genome_type,
                          docker_image=None):
    '''
    Run script to calculate correct reads

    :param tumour_wig: temporary wig file for tumour from readCounter
    :param normal_wig: temporary wig file for normal from readCounter
    :param outfile: temporary output text file (correct_reads.txt)
    :param config: dictionary of parameters for the run
    '''

    script = 'correct_reads.R'
    if not target_list:
        target_list = 'NULL'

    cmd = [script, tumour_wig, normal_wig, gc_wig, map_wig,
           target_list, outfile, genome_type]

    pypeliner.commandline.execute(*cmd)

    data = pd.read_csv(outfile, sep='\t', converters={'chr': str})
    if data['logR'].isnull().all():
        raise ValueError('all none')


def run_titan(
        infile, cnfile, outfile, obj_outfile, outparam,
        num_clusters, ploidy, sample_id, map_wig, titan_params, genome_type,
        docker_image=None, threads=8
):
    script = 'titan.R'

    cmd = [script, sample_id, infile, cnfile, map_wig, num_clusters,
           threads, ploidy, outfile, outparam,
           titan_params['myskew'], titan_params['estimate_ploidy'], titan_params['normal_param_nzero'],
           titan_params['normal_estimate_method'], titan_params['max_iters'], titan_params['pseudo_counts'],
           titan_params['txn_exp_len'], titan_params['txn_z_strength'], titan_params['alpha_k'],
           titan_params['alpha_high'], titan_params['max_copynumber'],
           titan_params['symmetric'], obj_outfile, genome_type, titan_params['chrom'],
           titan_params['y_threshold'], titan_params['max_depth']]

    pypeliner.commandline.execute(*cmd)


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def plot_titan(obj_file, output, tempdir, num_clusters, ploidy, chromosomes=None):
    if chromosomes is None:
        chromosomes = map(str, range(1, 23)) + ['X']

    script = 'plot_titan.R'

    chrom = 'c(' + ','.join(['\'' + c + '\'' for c in chromosomes]) + ')'

    try:
        cmd = [script, obj_file, tempdir, num_clusters, chrom, ploidy]
        pypeliner.commandline.execute(*cmd)

        cluster_ploidy_tempdir = os.path.join(tempdir, 'cluster_{}_ploidy_{}'.format(num_clusters, ploidy))
        pdfutils.merge_titan_pngs(cluster_ploidy_tempdir, output, num_clusters, chromosomes)
    except FileNotFoundError:
        chrom = '"{}"'.format(chrom)
        cmd = [script, obj_file, tempdir, num_clusters, chrom, ploidy]
        pypeliner.commandline.execute(*cmd)

        cluster_ploidy_tempdir = os.path.join(tempdir, 'cluster_{}_ploidy_{}'.format(num_clusters, ploidy))
        pdfutils.merge_titan_pngs(cluster_ploidy_tempdir, output, num_clusters, chromosomes)


def calc_cnsegments_titan(infile, outigv, outfile, sample_id):
    script = 'createTITANsegmentfiles.pl'

    symmetric = '1'

    cmd = [script, '-id=' + sample_id, '-infile=' + infile,
           '-outfile=' + outfile, '-outIGV=' + outigv,
           '-symmetric=' + symmetric]

    pypeliner.commandline.execute(*cmd)


def annot_pygenes(infile, outfile, pygenes_gtf):
    annotator = PygeneAnnotation(infile, outfile, gtf=pygenes_gtf)
    annotator.write_output()


def merge_to_h5(inputs, output, intervals, dtype=None):
    """

    :param inputs:
    :type inputs:
    :param output:
    :type output:
    :param intervals:
    :type intervals:
    :return:
    :rtype:
    """

    with pd.HDFStore(output, 'w', complevel=9, complib='blosc') as h5output:
        for interval in intervals:
            num_clusters = interval['num_clusters']
            ploidy = interval['ploidy']
            input_file = inputs[(num_clusters, ploidy)]

            input_df = pd.read_csv(input_file, sep='\t', dtype=dtype)

            tablename = '/cluster_{}/ploidy_{}'.format(num_clusters, ploidy)
            h5output.put(tablename, input_df, format='table')


def parse_titan_data(infile, titan_file, output):
    """
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    """

    parse_titan.parser(titan_file, infile, output)


def get_param_file_vals(params_file):
    '''
    exracts the values from a titan params
    file and returns them as a dict.

    :param params_file: titan params file
    '''
    params_file_vals = {}

    with open(params_file, 'r') as f:
        for line in f:
            k, v = line.split(":")
            params_file_vals[k] = np.array(v.split()).astype(float)

    return params_file_vals


def select_optimal_solution(
        chunks,
        params_files,
        segments,
        igv_segs,
        markers,
        parsed_files,
        plots,
        optimal_segment,
        optimal_igv_segs,
        optimal_param,
        optimal_marker,
        optimal_parsed,
        optimal_plot,
):
    '''
    selects the optimal cluster and ploidy
    combination from an input set of cluster/ploidy-
    resolved params and writes the corresponding
    set of optimal output files to an 'optimal'
    output directory.

    :params nclusts: number of clusters / sample
    :params nploidy: ploidy options/cluster/sample
    :params params_files: set of paramater files
        per ploidy/cluster
    :params segments: input set of output segments
        files per ploidy/cluster
    :params params: input set of output params
        files per ploidy/cluster
    :params markers: input set of output markers
        files per ploidy/cluster
    :params parsed_files: input set of parsed
        files per ploidy/cluster
    :params plots: input set of plots
        files per ploidy/cluster
    :params optimal_segment: output path
        for the optimal segment file
    :params optimal_param: output path
        for the optimal param file
    :params optimal_marker: output path
        for optimal marker file
    :params optimal_parsed: output path
        for optimal parsed file
    :params optimal_plots: output path
        for optimal plots file
    '''

    model_select_idxs = []

    # find optimal cluster/ploidy
    for chunk in chunks:
        params = params_files[chunk]
        parsed_params = get_param_file_vals(params)
        dbw_index = parsed_params['S_Dbw validity index (Both)'][0]
        model_select_idxs.append((chunk, dbw_index))

    best_model = min(model_select_idxs, key=lambda t: t[1])
    best_model = best_model[0]

    # copy the file at the optimal cluster/ploidy to the
    # optimal file output path
    csvutils.finalize_csv(segments[best_model], optimal_segment, sep='\t')
    csvutils.finalize_csv(params_files[best_model], optimal_param, sep='\t')
    csvutils.finalize_csv(markers[best_model], optimal_marker, sep='\t')
    csvutils.finalize_csv(parsed_files[best_model], optimal_parsed, sep='\t')
    shutil.copyfile(plots[best_model], optimal_plot)
    shutil.copyfile(igv_segs[best_model], optimal_igv_segs)

    with helpers.GetFileHandle(optimal_param, 'at') as params_output:
        ploidy, num_clusters = best_model
        params_output.write('ploidy: {}\n'.format(ploidy))
        params_output.write('num_clusters: {}\n'.format(num_clusters))


def tar_all_data(params, segs, igv_segs, markers, parsed, plots, tar_output, tempdir, chunks):
    helpers.makedirs(tempdir)

    for chunk in chunks:
        num_cluster, ploidy = chunk

        num_cluster = str(num_cluster)
        ploidy = str(ploidy)

        outdir = os.path.join(tempdir, 'numcluster_' + num_cluster, 'ploidy_' + ploidy)

        helpers.makedirs(outdir)

        params_outfile = os.path.join(outdir, 'params.csv')
        shutil.copyfile(params[chunk], params_outfile)

        segs_outfile = os.path.join(outdir, 'segs.csv')
        shutil.copyfile(segs[chunk], segs_outfile)

        igv_segs_outfile = os.path.join(outdir, 'igv_segs.csv')
        shutil.copyfile(igv_segs[chunk], igv_segs_outfile)

        markers_outfile = os.path.join(outdir, 'titan_markers.csv')
        shutil.copyfile(markers[chunk], markers_outfile)

        parsed_outfile = os.path.join(outdir, 'parsed.csv')
        shutil.copyfile(parsed[chunk], parsed_outfile)

        plots_outfile = os.path.join(outdir, 'plots.pdf')
        shutil.copyfile(plots[chunk], plots_outfile)

    helpers.make_tarfile(tar_output, tempdir)
