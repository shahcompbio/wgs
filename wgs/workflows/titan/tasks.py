import os
import tarfile

import pandas as pd
import pypeliner
import pysam
from wgs.utils import pdfutils
from wgs.utils import vcfutils

from scripts import PygeneAnnotation
from scripts import ReadCounter
from scripts import TransformVcfCounts


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            intervals.append(name + "_" + str(int(i * size)) + "_" + str(int((i + 1) * size)))

    return intervals


def merge_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)


def convert_museq_vcf2counts(infile, outfile, config):
    """
    Transform museq vcf file to a text file of counts

    :param infile: merged temporary vcf file from museq run (museq.vcf)
    :param outfile: temporary text file of counts (museq_postprocess.txt)
    """

    transformer = TransformVcfCounts(infile, outfile, config['dbsnp_positions'])
    transformer.main()


def run_readcounter(input_bam, output_wig, config):
    rc = ReadCounter(
        input_bam, output_wig, config['readcounter']['w'],
        config['chromosomes'], config['readcounter']['q'],
    )
    rc.main()


def calc_correctreads_wig(tumour_wig, normal_wig, target_list, outfile, config, docker_image=None):
    '''
    Run script to calculate correct reads

    :param tumour_wig: temporary wig file for tumour from readCounter
    :param normal_wig: temporary wig file for normal from readCounter
    :param outfile: temporary output text file (correct_reads.txt)
    :param config: dictionary of parameters for the run
    '''

    script = 'correctReads.R'
    gc = config['correction']['gc']
    map_wig = config['titan_params']['map']
    if not target_list:
        target_list = 'NULL'
    genome_type = config['titan_params']['genome_type']

    cmd = [script, tumour_wig, normal_wig, gc, map_wig,
           target_list, outfile, genome_type]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_titan(
        infile, cnfile, outfile, obj_outfile, outparam,
        titan_params, num_clusters, ploidy, sample_id,
        docker_image=None
):
    script = 'titan.R'
    map_wig = titan_params['map']

    cmd = [script, sample_id, infile, cnfile, map_wig, num_clusters,
           titan_params['num_cores'], ploidy, outfile, outparam,
           titan_params['myskew'], titan_params['estimate_ploidy'], titan_params['normal_param_nzero'],
           titan_params['normal_estimate_method'], titan_params['max_iters'], titan_params['pseudo_counts'],
           titan_params['txn_exp_len'], titan_params['txn_z_strength'], titan_params['alpha_k'],
           titan_params['alpha_high'], titan_params['max_copynumber'],
           titan_params['symmetric'], obj_outfile, titan_params['genome_type'], titan_params['chrom'],
           titan_params['y_threshold'], titan_params['max_depth']]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def plot_titan(obj_file, output, tempdir, num_clusters, ploidy, docker_image=None):
    script = 'plot_titan.R'

    chrom = 'c(1:22,\'X\')'

    cmd = [script, obj_file, tempdir, num_clusters, chrom, ploidy]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    chromosomes = map(str, range(1, 23)) + ['X']
    pdfutils.merge_pngs(tempdir, output, chromosomes)


def calc_cnsegments_titan(infile, outigv, outfile, docker_image=None):
    script = 'createTITANsegmentfiles.pl'

    sample_id = os.path.basename(infile).split('_')[0]
    symmetric = '1'

    cmd = [script, '-id=' + sample_id, '-infile=' + infile,
           '-outfile=' + outfile, '-outIGV=' + outigv,
           '-symmetric=' + symmetric]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def annot_pygenes(infile, outfile, config):
    gene_sets_gtf = config['pygenes_gtf']
    annotator = PygeneAnnotation(infile, outfile, gtf=gene_sets_gtf)
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


def parse_titan(infile, params_file, titan_file, output, config, sample_id, docker_image=None):
    """
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    """

    cmd = [
        'vizutils_parse_titan',
        '--infile', infile,
        '--paramsfile', params_file,
        '--titanfile', titan_file,
        '--output', output,
        '--case_id', sample_id,
        '--tumour_id', sample_id,
        '--normal_id', sample_id + 'N',
    ]

    for key, val in config.iteritems():
        if val is None:
            continue
        elif isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
