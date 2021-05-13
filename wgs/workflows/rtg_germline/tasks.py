'''
Created on Feb 21, 2018

@author: dgrewal
'''
import os
import shutil

import pypeliner
import pysam
from wgs.utils import bamutils
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


def rtg_germline_command(reference, interval, bam_file, tempdir):
    interval = interval.split('_')

    interval = '{}:{}-{}'.format(interval[0], interval[1], interval[2])

    cmd = [
        'rtg', 'snp', '-o', tempdir, '-t', reference, bam_file,
        '--no-calibration', '-m', 'illumina',
        '--region', interval
    ]

    return cmd


def rtg_move_vcf(tempdir, vcf_output):
    if '.vcf.gz' not in vcf_output:
        raise Exception('RTG vcf files are gzipped')

    vcf_file = os.path.join(tempdir, 'snps.vcf.gz')
    shutil.copyfile(vcf_file, vcf_output)

    tbi_file = os.path.join(tempdir, 'snps.vcf.gz.tbi')
    shutil.copyfile(tbi_file, vcf_output + '.tbi')


def run_rtg_germline(vcf, reference, interval, bam_file, tempdir):
    # rtg fails if output dir already exists
    helpers.rmdirs(tempdir)

    cmd = rtg_germline_command(reference, interval, bam_file, tempdir)
    pypeliner.commandline.execute(*cmd)

    vcf_file = os.path.join(tempdir, 'snps.vcf.gz')

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(vcf_file, vcf, normal_id)


def run_rtg_one_job(
        tempdir, vcf, reference, intervals, bam_file, freebayes_docker_image=None, vcftools_docker_image=None
):
    helpers.rmdirs(tempdir)
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        cmd = rtg_germline_command(reference, interval, bam_file, ival_temp_dir)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, freebayes_docker_image)

    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        output = os.path.join(ival_temp_dir, 'germline.vcf.gz')
        rtg_move_vcf(ival_temp_dir, output)

    vcf_files = [os.path.join(tempdir, str(i), 'germline.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'germline_merge')
    helpers.makedirs(merge_tempdir)

    temp_vcf = os.path.join(merge_tempdir, 'merged_rtg.vcf')
    merge_vcfs(vcf_files, temp_vcf, merge_tempdir, docker_image=vcftools_docker_image)

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(temp_vcf, vcf, normal_id)


def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)
