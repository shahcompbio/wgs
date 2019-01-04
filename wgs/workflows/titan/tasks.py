import os
import pypeliner
import multiprocessing
import subprocess
from wgs.utils import vcfutils, helpers
import pysam
import tarfile
import pandas as pd

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length/size)+1)):
            intervals.append( name+ "_" + str(i*size) +"_"+ str((i+1)*size))

    return intervals

def run_museq(tumour_bam, normal_bam, out, log, reference, interval, museq_params):
    '''
    Run museq script for all chromosomes and merge vcf files

    :param tumour_bam: path to tumour bam
    :param tumour_bai: path to tumour bai
    :param normal_bam: path to normal bam
    :param normal_bai: path to normal bai
    :param out: temporary output vcf file for the merged vcf files (museq.vcf)
    :param log: temporary log file (museq.log)
    :param config: dictionary of parameters for the run
    '''
    interval = interval.split('_')
    interval = interval[0] + ':' + interval[1] + '-' + interval[2]

    cmd = ['museq_het','normal:' + normal_bam, 'tumour:' + tumour_bam,
           'reference:' + reference, '--out', out,
           '--log', log, '--interval', interval, '--verbose']

    for key, val in museq_params.iteritems():
        if isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd)



def merge_vcfs(inputs, output):
    vcfutils.concatenate_vcf(inputs, output)

def convert_museq_vcf2counts(infile, outfile, config):
    '''
    Transform museq vcf file to a text file of counts

    :param infile: merged temporary vcf file from museq run (museq.vcf)
    :param outfile: temporary text file of counts (museq_postprocess.txt)
    :param config: dictionary of parameters for the run
    '''

    script = os.path.join(scripts_directory, 'transform_vcf_to_counts.py')
    positions_file = config['dbsnp_positions']

    cmd = ['python', script, '--infile', infile, '--outfile', outfile,
            '--positions_file', positions_file]

    pypeliner.commandline.execute(*cmd)


def run_readcounter(bam, outfile, config):
    '''
    Run readCounter on the input bam file

    :param bam: path to normal or tumour bam
    :param bai: path to normal or tumour bai
    :param outfile: temporary wig file
    :param config: dictionary of parameters for the run
    '''

    readcounter = 'readCounter'
    w = config['readcounter']['w']
    q = config['readcounter']['q']
    chrs = ','.join(config['chromosomes'])

    cmd = [readcounter, '-w', w, '-q', q, '-c', chrs, bam, '>', outfile]

    pypeliner.commandline.execute(*cmd)


def calc_correctreads_wig(tumour_wig, normal_wig, outfile, config):
    '''
    Run script to calculate correct reads

    :param tumour_wig: temporary wig file for tumour from readCounter
    :param normal_wig: temporary wig file for normal from readCounter
    :param outfile: temporary output text file (correct_reads.txt)
    :param config: dictionary of parameters for the run
    '''

    script = os.path.join(scripts_directory, 'correctReads.R')
    gc = config['correction']['gc']
    map_wig = config['titan_params']['map']
    target_list = 'NULL'
    genome_type = config['titan_params']['genome_type']

    cmd = ['Rscript', script, tumour_wig, normal_wig, gc, map_wig, 
            target_list, outfile, genome_type]

    pypeliner.commandline.execute(*cmd)


def run_titan(infile, cnfile, outfile, obj_outfile, outparam, titan_params, num_clusters, ploidy):
    script = os.path.join(scripts_directory, 'titan.R')
    map_wig = titan_params['map']

    sample_id = os.path.basename(outfile).split('_')[0]

    cmd = ['Rscript', script, sample_id, infile, cnfile, map_wig, num_clusters,
           titan_params['num_cores'], ploidy, outfile, outparam,
           titan_params['myskew'], titan_params['estimate_ploidy'], titan_params['normal_param_nzero'],
           titan_params['normal_estimate_method'], titan_params['max_iters'], titan_params['pseudo_counts'],
           titan_params['txn_exp_len'], titan_params['txn_z_strength'], titan_params['alpha_k'],
           titan_params['alpha_high'], titan_params['max_copynumber'],
           titan_params['symmetric'], obj_outfile, titan_params['genome_type'], titan_params['chrom'],
           titan_params['y_threshold'], titan_params['max_depth']]

    pypeliner.commandline.execute(*cmd)

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

def plot_titan(obj_file, titan_input, output, tempdir, config, num_clusters, ploidy):
    script = os.path.join(scripts_directory, 'plot_titan.R')

    chrom = 'c(1:22,\'X\')'

    cmd = ['Rscript', script, obj_file, tempdir, num_clusters, chrom, ploidy]

    pypeliner.commandline.execute(*cmd)

    make_tarfile(output, tempdir)

def calc_cnsegments_titan(infile, outigv, outfile):
    script = os.path.join(scripts_directory, 'createTITANsegmentfiles.pl')

    sample_id = os.path.basename(infile).split('_')[0]
    symmetric = '1'

    cmd = ['perl', script, '-id=' + sample_id, '-infile=' + infile,
            '-outfile=' + outfile, '-outIGV=' + outigv, 
            '-symmetric=' + symmetric]

    helpers.run_cmd(cmd)


def annot_pygenes(infile, outfile, config):
    script = os.path.join(scripts_directory, 'pygene_annotation.py')
    gene_sets_gtf = config['pygenes_gtf']

    cmd = ['python', script, '--infile', infile, '--outfile', outfile,
            '--gene_sets_gtf', gene_sets_gtf]

    helpers.run_cmd(cmd)



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