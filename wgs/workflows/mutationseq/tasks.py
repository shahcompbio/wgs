'''
Created on Feb 21, 2018

@author: pwalters
'''
import os
import pysam
import pypeliner
import multiprocessing
from wgs.utils import helpers, vcfutils


scripts_directory = os.path.join(
    os.path.realpath(os.path.dirname(__file__)), 'scripts')


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

def run_museq(out, log, reference, interval, museq_params, tumour_bam=None,
              normal_bam=None):
    '''
    Run museq script for all chromosomes and merge VCF files

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to the temporary output VCF file for the merged VCF files
    :param log: path to the log file
    :param config: path to the config YAML file
    '''

    cmd = ['museq']

    if tumour_bam:
        cmd.append('tumour:' + tumour_bam)
    if normal_bam:
        cmd.append('normal:' + normal_bam)

    interval = interval.split('_')
    interval = interval[0] + ':' + interval[1] + '-' + interval[2]

    cmd.extend(['reference:' + reference, '--out', out,
                '--log', log, '--interval', interval, '-v'])

    if not tumour_bam or not normal_bam:
        cmd.extend(['-s'])

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

def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)

def run_museqportrait(infile, out_pdf, out_txt, museqportrait_log, single_mode, config, sample):
    '''
    Run museqportrait script on the input VCF file

    :param infile: temporary input VCF file
    :param out_dir: temporary output VCF file
    :param museqportrait_log: path to the log file
    '''

    if single_mode:
        script = os.path.join(scripts_directory, 'portraits_single_sample.py')
        cmd = ['python', script, '--thousand_gen', config["thousandgen_params"]["db"],
               '--output', out_pdf, '--data', out_txt, '--threshold', config['threshold'],
               '--dbsnp', config["dbsnp_params"]["db"], '--log_file', museqportrait_log,
               '--ref_data', config['refdata_single_sample'], '--variant_label',
               sample, '--variant_file', infile, '--tabix_path', 'tabix'
        ]
        pypeliner.commandline.execute(*cmd)

    else:
        cmd = ['museqportrait', '--log', museqportrait_log, '--output-pdf',
               out_pdf, '--output-txt', out_txt, infile]
        pypeliner.commandline.execute(*cmd)
