'''
Created on Feb 27, 2018

@author: dgrewal
'''
import os
import warnings

import pypeliner
from wgs.utils import helpers

def _get_header(infile):
    '''
    Extract header from the VCF file

    :param infile: input VCF file
    :return: header
    '''

    header = []
    for line in infile:
        if line.startswith('##'):
            header.append(line)
        elif line.startswith('#'):
            header.append(line)
            return header
        else:
            raise Exception('invalid header: missing #CHROM line')

    warnings.warn("One of the input files is empty")
    return []


def concatenate_vcf(infiles, outfile):
    '''
    Concatenate VCF files

    :param infiles: dictionary of input VCF files to be concatenated
    :param outfile: output VCF file
    '''
    if isinstance(infiles, dict):
        keys = infiles.keys()
        keys = sorted(keys)
        infiles = [infiles[val] for val in keys]

    with helpers.GetFileHandle(outfile, 'w') as ofile:
        header = None

        for ifile in infiles:

            if os.path.getsize(ifile) == 0:
                warnings.warn('input file {} is empty'.format(ifile))
                continue

            with helpers.GetFileHandle(ifile) as f:

                if not header:
                    header = _get_header(f)

                    for line in header:
                        ofile.write(line)
                else:
                    if not _get_header(f) == header:
                        warnings.warn('merging vcf files with mismatching headers')

                for l in f:
                    ofile.write(l)


def sort_vcf(infile, outfile, docker_image=None):
    cmd = ['cat', infile, '|', 'vcf-sort', '>', outfile]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def update_germline_header_sample_ids(infile, outfile, sample_id):
    with helpers.GetFileHandle(infile) as indata:
        with helpers.GetFileHandle(outfile, 'wt') as outdata:
            for line in indata:
                if line.startswith('#CHROM'):
                    outdata.write('##normal_sample={}\n'.format(sample_id))
                    line = line.replace('NORMA', sample_id)
                    outdata.write(line)
                else:
                    outdata.write(line)
