'''
Created on Feb 21, 2018

@author: pwalters
'''
import os

import pypeliner
import pysam
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
        for i in range(int((length / size) + 1)):
            start = str(int(i * size))
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def run_museq(out, log, reference, interval, museq_params, tumour_bam=None,
              normal_bam=None, return_cmd=False):
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
    if len(interval) == 1:
        interval = interval[0]
    elif len(interval) == 2:
        interval = interval[0] + ':' + interval[1]
    elif len(interval) == 3:
        interval = interval[0] + ':' + interval[1] + '-' + interval[2]
    else:
        raise Exception()

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

    if return_cmd:
        return cmd
    else:
        pypeliner.commandline.execute(*cmd)


def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)


def run_museqportrait(infile, out_pdf, out_txt, museqportrait_log, single_mode, config):
    '''
    Run museqportrait script on the input VCF file

    :param infile: temporary input VCF file
    :param out_dir: temporary output VCF file
    :param museqportrait_log: path to the log file
    '''

    if single_mode:
        script = os.path.join(scripts_directory, 'singlesampleplot.py')
        cmd = ['python', script, '--thousand_gen', config["thousandgen_params"]["db"],
               '--output', out_pdf, '--threshold', config['threshold'],
               '--dbsnp', config["dbsnp_params"]["db"],
               '--ref_data', config['refdata_single_sample'],
               '--variant_file', infile,
               ]
        pypeliner.commandline.execute(*cmd)

        # touch the txt file to avoid pypeliner errors
        open(out_txt, 'w').close()
        open(museqportrait_log, 'w').close()

    else:
        cmd = ['museqportrait', '--log', museqportrait_log, '--output-pdf',
               out_pdf, '--output-txt', out_txt, infile]
        pypeliner.commandline.execute(*cmd)




def run_museq_one_job(tempdir, museq_vcf, reference, intervals, museq_params, tumour_bam=None,
              normal_bam=None):
    '''
    Run museq script for all chromosomes and merge VCF files

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to the temporary output VCF file for the merged VCF files
    :param log: path to the log file
    :param config: path to the config YAML file
    '''

    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output = os.path.join(ival_temp_dir, 'museq.vcf')
        log = os.path.join(ival_temp_dir, 'museq.log')

        command = run_museq(
            output, log, reference, interval, museq_params,
            tumour_bam=tumour_bam, normal_bam=normal_bam,
            return_cmd=True
        )

        commands.append(command)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, None)

    vcf_files = [os.path.join(tempdir, str(i), 'museq.vcf') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'museq_merge')
    helpers.makedirs(merge_tempdir)
    merge_vcfs(vcf_files, museq_vcf, merge_tempdir)