'''
Created on 30 Apr 2015

@author: jrosner
@last_update: 30 Apr 2015
'''

import os
import warnings

version = '1.0.1'


def get_table_hdr():
    '''
    return the header from a
    mutation assessor file as a list
    '''

    ma_file = os.path.join(args.db, 'MA.chr01.txt')
    with open(ma_file) as f:
        line = f.readline()

    return line.rstrip().split('\t')


def load_table(chrom):
    ''' load mutationassessor data for chromosome chrom '''

    print 'loading table for chromosome {}'.format(chrom)
    ma_table = {}
    chrom = chrom.replace('chr','')

    if chrom != 'X' and chrom != 'Y':
        chrom = chrom.zfill(2)

    ma_file = os.path.join(args.db, 'MA.chr{}.txt'.format(chrom))
    with open(ma_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip().split()
            key = '_'.join(line[0].split(',')[1:5])
            ma_table[key] = line

    return ma_table


def get_ma_description():
    ''' add new mutation assessor header line '''

    ma_vcf_hdr = ('##INFO=<ID=MA,Number=.,Type=String,Description='
                  '\"Predicted functional impact of amino-acid substitutions'
                  'in proteins.Format: ({}) \">\n'.format('|'.join(get_table_hdr())))
    return ma_vcf_hdr


def write_hdr(infile, o):
    ''' write the header to the output file '''

    ma_line = get_ma_description()
    hdr = []
    with open(infile, 'r') as i:
        for line in i:
            if line.startswith('#'):
                hdr.append(line)
            else:
                break

    if hdr:
        hdr = hdr[:-1] + [ma_line] + [hdr[-1]]
        o.write(''.join(hdr))


def annot_lookup(key, table):
    ''' lookup and format annotation '''

    annot = table.get(key)
    if annot:
        print 'found match: {}'.format(key)
        annot = ';MA=({})'.format('|'.join(annot))
    else:
        annot = ';MA=()'

    return annot


def annot_insert(l, a):
    ''' insert ma annotation into the info section of the vcf line '''

    # info section in vcf is column 7 (zero-based)
    info    = 7
    l[info] = l[info] + a

    return '\t'.join(l)


def write_annot_data(infile, out):
    ''' write the annotated data to out '''

    ma_table     = None
    chr_seen     = []
    chr_no_table = []
    chr_current  = None
    chr_table    = None

    with open(infile, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue

            line        = line.rstrip().split('\t')
            chr_current = line[0]

            if chr_current in chr_no_table:
                continue

            if chr_table is None or chr_current != chr_table:
                if chr_current in chr_seen:
                    warnings.warn('detected unsorted bam file, annotation could be very slow')
                    print chr_current
                    print chr_table
                    print chr_seen

                try:
                    ma_table  = load_table(chr_current)
                except IOError:
                    print 'no matching table found for {}'.format(chr_current)
                    chr_no_table.append(chr_current)
                    out.write('\t'.join(line) + '\n')
                    continue

                chr_table = chr_current
                chr_seen.append(chr_current)

            pos   = line[1]
            ref   = line[3]
            alt   = line[4]
            key   = '_'.join([chr_current, pos, ref, alt])
            annot = annot_lookup(key, ma_table)
            line  = annot_insert(line, annot) + '\n'

            out.write(line)


def _main():
    ''' main function '''

    with open(args.output,'w') as out:

        write_hdr(args.vcf, out)
        write_annot_data(args.vcf, out)



if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf',
                        required=True,
                        help='input vcf file')
    parser.add_argument('--output',
                        required = True,
                        help = 'output vcf file')
    parser.add_argument('--db',
                        required = True,
                        help = 'mutation assessor db directory')

    args = parser.parse_args()

    _main()

