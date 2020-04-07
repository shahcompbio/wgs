#!/usr/bin/env python

import argparse
import gzip
import os
import sys
# from string import maketrans

import pysam
from wgs.utils import helpers


class OpenFile(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.filename.endswith(".gz"):
            self.handle = gzip.open(self.filename, self.mode)
        else:
            self.handle = open(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()


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


def get_bam_reader(infile):
    if infile.endswith('.sam'):
        bam = pysam.Samfile(infile, 'r', check_sq=False)
    else:
        bam = pysam.Samfile(infile, "rb", check_sq=False)

    return bam


def write_to_fastq(aln, read, rg, out_fastq):
    if aln.is_reverse:
        outstr = "@" + str(aln.qname) + "/" + str(read) + " " + "RG:Z:" + str(rg) + "\n" + str(
            revcomp(aln.seq)) + "\n" + "+" + "\n" + str(aln.qual[::-1])
    else:
        outstr = "@" + str(aln.qname) + "/" + str(read) + " " + "RG:Z:" + str(rg) + "\n" + str(
            aln.seq) + "\n" + "+" + "\n" + str(aln.qual)

    out_fastq.write(outstr + '\n')


def revcomp(seq):
    seq1 = seq.translate(str.maketrans("AGCTagct", "TCGAtcga"))
    seq2 = seq1[::-1]
    return seq2


def get_outfiles(outdir, readgroups):
    outfiles = {}

    for readgroup in readgroups:
        helpers.makedirs(os.path.join(outdir, readgroup))
        r1 = os.path.join(outdir, readgroup, 'R1.fastq.gz')
        r2 = os.path.join(outdir, readgroup, 'R2.fastq.gz')
        outfiles[readgroup] = (r1, r2)

    return outfiles


def open_outfiles(outfiles):
    opened_files = {}
    for rgid, fastqs in outfiles.items():
        opened_files[rgid] = (
            OpenFile(fastqs[0], 'wt').__enter__(),
            OpenFile(fastqs[1], 'wt').__enter__()
        )
    return opened_files


def close_outfiles(outfiles):
    for rgid, fastqs in outfiles.items():
        fastqs[0].close()
        fastqs[1].close()


def bam_to_fastq(infile, outdir):
    readgroups = get_read_groups(infile)

    outfiles = get_outfiles(outdir,readgroups)
    outfiles = open_outfiles(outfiles)

    read_data = {}

    bam = get_bam_reader(infile)

    for al in bam:

        # must be primary read alignment, not secondary or supplementary
        if al.is_secondary or al.flag & 2048 == 2048:
            continue

        # skip unpaired reads
        if not al.is_paired:
            continue

        # ensures the read is not hard-clipped. important
        # when the BAM doesn't have shorter hits flagged as
        # secondary
        if al.cigar is not None and 5 in [x[0] for x in al.cigar]:
            continue

        # add read name to dictionary if not already there
        key = al.qname
        if key not in read_data:
            read_data.setdefault(key, al)
        # print matched read pairs
        else:
            # RG:Z:ID
            try:
                RG1 = read_data[key].opt('RG')
            except KeyError:
                RG1 = ""
            try:
                RG2 = al.opt('RG')
            except KeyError:
                RG2 = ""

            rgid = str(al.opt('RG'))
            fastq_r1 = outfiles[rgid][0]
            fastq_r2 = outfiles[rgid][1]

            if al.is_read1:
                write_to_fastq(al, 1, RG2, fastq_r1)
                write_to_fastq(read_data[key], 2, RG1, fastq_r2)
            else:
                write_to_fastq(read_data[key], 1, RG1, fastq_r1)
                write_to_fastq(al, 2, RG2, fastq_r2)
            del read_data[key]
    if len(read_data) != 0:
        sys.stderr.write('Warning: %s unmatched name groups\n' % len(read_data))

    close_outfiles(outfiles)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input',
        help='input bam file'
    )

    parser.add_argument(
        'outdir',
        help='output directory'
    )

    args = parser.parse_args()

    args = vars(args)
    return args


def main():
    args = parse_args()
    bam_to_fastq(args['input'], args['outdir'])


if __name__ == '__main__':
    main()
