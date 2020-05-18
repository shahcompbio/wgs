'''
Created on Feb 21, 2018

@author: pwalters
'''
import os

import pypeliner
import pysam
from wgs.utils import helpers
from wgs.utils import vcfutils


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


def samtools_germline_command(vcf, reference, interval, bam_file):
    interval = interval.split('_')

    interval = '{}:{}-{}'.format(interval[0], interval[1], interval[2])

    cmd = [
        'samtools', 'mpileup', '-ugf', reference, '-Q', 20, '-q', 10, '-r', interval, bam_file,
        '|', 'bcftools', 'call', '-vmO', 'z', '-o', vcf
    ]

    return cmd


def run_samtools_germline(vcf, reference, interval, bam_file, docker_image=None):
    cmd = samtools_germline_command(vcf, reference, interval, bam_file)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_samtools_germline_one_job(
        tempdir, vcf, reference, intervals, bam_file, samtools_docker_image=None, vcftools_docker_image=None
):
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output = os.path.join(ival_temp_dir, 'germline.vcf.gz')
        cmd = samtools_germline_command(output, reference, interval, bam_file)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, samtools_docker_image)

    vcf_files = [os.path.join(tempdir, str(i), 'germline.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'germline_merge')
    helpers.makedirs(merge_tempdir)
    merge_vcfs(vcf_files, vcf, merge_tempdir, docker_image=vcftools_docker_image)


def merge_vcfs(inputs, outfile, tempdir, docker_image=None):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile, docker_image=docker_image)


def roh_calling(samtools_germlines, roh_output, docker_image=None):
    cmd = ['bcftools', 'roh', '-G30', '--AF-dflt', 0.4, samtools_germlines, '>', roh_output]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
