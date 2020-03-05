'''
Created on Aug, 2019

@author: Diljot Grewal
'''

import os

import pypeliner

from wgs.utils import helpers

import pysam


def split_bam_by_readgroups(infile, readgroups, r1_fastq_output, r2_fastq_output):
    for readgroup in readgroups:
        cmd = ['wgs_bamtofastq', infile, r1_fastq_output[readgroup], r2_fastq_output[readgroup], '--readgroup',
               readgroup]
        pypeliner.commandline.execute(*cmd)


def split_by_rg(infile, read1_output, read2_output, tempdir):

    helpers.makedirs(tempdir)

    cmd = ['wgs_bamtofastq', infile, tempdir]
    pypeliner.commandline.execute(*cmd)

    readgroups = os.listdir(tempdir)

    for readgroup in readgroups:

        os.rename(
            os.path.join(tempdir, readgroup, 'R1.fastq.gz'),
            read1_output[readgroup]
        )

        os.rename(
            os.path.join(tempdir, readgroup, 'R2.fastq.gz'),
            read2_output[readgroup]
        )


def get_read_group(infile):
    bam = pysam.AlignmentFile(infile, mode='rb', check_sq=False)
    header = bam.header['RG']
    return {rg['ID']: rg for rg in header}
