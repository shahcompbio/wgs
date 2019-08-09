
'''
Created on Aug, 2019

@author: Diljot Grewal
'''

import pypeliner
import pysam


def get_read_groups(file_name):
    config = dict()

    bam = pysam.AlignmentFile(file_name, mode='rb', check_sq=False)

    readgroup_ref = {info['ID']: info for info in bam.header['RG']}

    header_text = bam.text

    for line in header_text.split('\n'):
        if line.startswith("@RG"):
            rgid = line.split()[1]
            assert rgid.startswith("ID:")
            rgid = rgid[3:]
            config[rgid] = (readgroup_ref[rgid], line)

    assert len(config) == len(readgroup_ref)
    bam.close()
    return config


def split_bam_by_readgroups(infile, readgroups, r1_fastq_output, r2_fastq_output):
    for readgroup in readgroups:
        cmd = ['wgs_bamtofastq', infile,r1_fastq_output[readgroup], r2_fastq_output[readgroup] ,'--readgroup',  readgroup ]
        pypeliner.commandline.execute(*cmd)

def split_by_rg(infile, read1_output, read2_output):
    readgroups = get_read_groups(infile)

    split_bam_by_readgroups(infile, readgroups, read1_output, read2_output)

