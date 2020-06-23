'''
Created on Feb 21, 2018

@author: dgrewal
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
            end = str(min(int((i + 1) * size), length))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def get_sample_id(bamfile):
    bam = pysam.AlignmentFile(bamfile)
    readgroups = bam.header['RG']

    samples = set()

    for readgroup in readgroups:
        samples.add(readgroup['SM'])

    assert len(samples) == 1

    return list(samples)[0]


def mutect_run_command(reference, interval, normal_bam, tumour_bam, vcf_out):
    interval = interval.split('_')

    interval = '{}:{}-{}'.format(interval[0], interval[1], interval[2])

    normal_sample_id = get_sample_id(normal_bam)

    cmd = [
        'gatk', 'Mutect2', '--input', normal_bam, '--input', tumour_bam,
        '--normal', normal_sample_id, '-O', vcf_out, '--intervals', interval,
        '-R', reference
    ]

    return cmd


def mutect_filter_command(reference, vcf_in, vcf_out):
    cmd = [
        'gatk', 'FilterMutectCalls', '-R', reference, '-V', vcf_in, '-O', vcf_out
    ]

    return cmd


def run_mutect(vcf, reference, interval, normal_bam, tumour_bam, tempdir, docker_image=None):
    helpers.makedirs(tempdir)
    unfiltered_vcf = os.path.join(tempdir, 'temp.vcf')
    cmd = mutect_run_command(reference, interval, normal_bam, tumour_bam, unfiltered_vcf)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cmd = mutect_filter_command(reference, unfiltered_vcf, vcf)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_mutect_one_job(
        tempdir, vcf, reference, intervals, normal_bam, tumour_bam, freebayes_docker_image=None,
        vcftools_docker_image=None
):
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        unfiltered_output = os.path.join(ival_temp_dir, 'mutect.vcf.gz')
        cmd = mutect_run_command(reference, interval, normal_bam, tumour_bam, unfiltered_output)
        commands.append(cmd)

        output = os.path.join(ival_temp_dir, 'mutect.vcf.gz')
        cmd = mutect_filter_command(reference, unfiltered_output, output)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, freebayes_docker_image)

    vcf_files = [os.path.join(tempdir, str(i), 'mutect.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'mutect_merge')
    helpers.makedirs(merge_tempdir)
    merge_vcfs(vcf_files, vcf, merge_tempdir, docker_image=vcftools_docker_image)


def merge_vcfs(inputs, outfile, tempdir, docker_image=None):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile, docker_image=docker_image)
