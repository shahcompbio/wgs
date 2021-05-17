'''
Created on Feb 21, 2018

@author: dgrewal
'''
import os

import pypeliner
import pysam
from wgs.utils import helpers
from wgs.utils import vcfutils
from wgs.utils import bamutils


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
            end = str(min(int((i + 1) * size), length))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def freebayes_germline_command(vcf, reference, interval, bam_file):
    interval = interval.split('_')

    interval = '{}:{}-{}'.format(interval[0], interval[1], interval[2])

    cmd = [
        'freebayes', '-f', reference, '-r', interval, bam_file, '>', vcf
    ]

    return cmd


def run_freebayes_germline(vcf, reference, interval, bam_file, tempdir):
    helpers.makedirs(tempdir)
    temp_vcf = os.path.join(tempdir, 'temp.vcf')

    cmd = freebayes_germline_command(temp_vcf, reference, interval, bam_file)
    pypeliner.commandline.execute(*cmd)

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(temp_vcf, vcf, normal_id)



def run_freebayes_one_job(
        tempdir, vcf, reference, intervals, bam_file,
):
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output = os.path.join(ival_temp_dir, 'germline.vcf.gz')
        cmd = freebayes_germline_command(output, reference, interval, bam_file)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir)

    vcf_files = [os.path.join(tempdir, str(i), 'germline.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'germline_merge')
    helpers.makedirs(merge_tempdir)
    temp_vcf = os.path.join(merge_tempdir,'temp_freebayes.vcf')
    merge_vcfs(vcf_files, temp_vcf, merge_tempdir)

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(temp_vcf, vcf, normal_id)


def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)
