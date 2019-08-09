#!/usr/bin/env python

import argparse
import sys
from string import *

import pysam
import gzip

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
    seq1 = seq.translate(maketrans("AGCTagct", "TCGAtcga"))
    seq2 = seq1[::-1]
    return seq2


def bam_to_fastq(infile, out_r1, out_r2, readgroup):
    read_data = {}

    bam = get_bam_reader(infile)

    with OpenFile(out_r1, 'w') as fastq_r1, OpenFile(out_r2, 'w') as fastq_r2:

        for al in bam:

            # must be primary read alignment, not secondary or supplementary
            if al.is_secondary or al.flag & 2048 == 2048:
                continue

            # skip unpaired reads
            if not al.is_paired:
                continue

            if not str(al.opt('RG')) == str(readgroup):
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

                if al.is_read1:
                    write_to_fastq(al, 1, RG2, fastq_r1)
                    write_to_fastq(read_data[key], 2, RG1, fastq_r2)
                else:
                    write_to_fastq(read_data[key], 1, RG1, fastq_r1)
                    write_to_fastq(al, 2, RG2, fastq_r2)
                del read_data[key]
        if len(read_data) != 0:
            sys.stderr.write('Warning: %s unmatched name groups\n' % len(read_data))


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input',
        help='input bam file'
    )

    parser.add_argument(
        'R1',
        help='input bam file'
    )

    parser.add_argument(
        'R2',
        help='input bam file'
    )

    parser.add_argument(
        '--readgroup',
        help='input bam file'
    )

    args = parser.parse_args()

    args = vars(args)
    return args


def main(args):
    args = parse_args()
    bam_to_fastq(args['input'], args['R1'], args['R2'], args['readgroup'])


if __name__ == '__main__':
    main(args)
