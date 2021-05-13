'''
Created on Feb 21, 2018

@author: dgrewal
'''
import os

import pypeliner
import pysam
from wgs.utils import bamutils
from wgs.utils import csvutils
from wgs.utils import helpers
from wgs.utils import vcfutils

import pandas as pd


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


def run_samtools_germline(vcf, reference, interval, bam_file, tempdir):
    helpers.makedirs(tempdir)
    vcf_file = os.path.join(tempdir, 'samtools_snps.vcf.gz')

    cmd = samtools_germline_command(vcf_file, reference, interval, bam_file)
    pypeliner.commandline.execute(*cmd)

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(vcf_file, vcf, normal_id)


def run_samtools_germline_one_job(
        tempdir, vcf, reference, intervals, bam_file
):
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output = os.path.join(ival_temp_dir, 'germline.vcf.gz')
        cmd = samtools_germline_command(output, reference, interval, bam_file)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir)

    vcf_files = [os.path.join(tempdir, str(i), 'germline.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'germline_merge')
    helpers.makedirs(merge_tempdir)

    temp_vcf = os.path.join(merge_tempdir, 'merged_rtg.vcf')
    merge_vcfs(vcf_files, temp_vcf, merge_tempdir)

    normal_id = bamutils.get_sample_id(bam_file)
    vcfutils.update_germline_header_sample_ids(temp_vcf, vcf, normal_id)


def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)


def roh_calling(samtools_germlines, roh_output, tempdir):
    helpers.makedirs(tempdir)

    output = os.path.join(tempdir, 'output.csv')

    cmd = ['bcftools', 'roh', '-G30', '--AF-dflt', 0.4, samtools_germlines, '>', output]

    pypeliner.commandline.execute(*cmd)

    parse_roh_output( output, roh_output)


def parse_roh_output(infile, outfile):
    parsed = []

    with helpers.GetFileHandle(infile) as indata:
        for line in indata:
            if line.startswith('#'):
                continue

            line = line.strip().split()

            if line[0] == 'ST':
               parsed.append({
                   'type': line[0],
                   'sample': line[1],
                   'chromosome': line[2],
                   'start': line[3],
                   'end': float('nan'),
                   'state': line[4],
                   'length': float('nan'),
                   'num_markers': float('nan'),
                   'quality': line[5]
               })
            elif line[0] == 'RG':
                parsed.append({
                    'type': line[0],
                    'sample': line[1],
                    'chromosome': line[2],
                    'start': line[3],
                    'end': line[4],
                    'state': float('nan'),
                    'length': line[5],
                    'num_markers': line[6],
                    'quality': line[7]
                })


    parsed = pd.DataFrame(parsed)

    csvutils.write_dataframe_to_csv_and_yaml(parsed, outfile, write_header=True)
