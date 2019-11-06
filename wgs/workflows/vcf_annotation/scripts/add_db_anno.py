'''
@author:  Celia Siu
created: 2014-07-31 (Thursday)
last_updated: 12 Feb 16 by dgrewal

'''
import argparse
import os
import warnings
from collections import defaultdict

import pysam

version = '1.3.1'


def resolve_db_position(input_type, pos):
    '''
    adjust database position to match convention of
    input file type for indexing a deletion
    '''
    if input_type == 'snv':
        pos += 1

    return pos


def load_chromosome(db, chromosome, input_type, flag_with_id):
    db_dict = defaultdict(list)

    f = pysam.Tabixfile(db)
    try:
        fetchdata = f.fetch(chromosome)
    except:
        warnings.warn("warning")
        return db_dict
    for line in fetchdata:

        chrom, pos = line.split('\t')[0:2]

        # resolve position when record is a deletion
        if len(line.split()[3]) > len(line.split()[4]):
            pos = resolve_db_position(input_type, int(pos))

        chrom.replace('chr', '')
        key = ','.join([chrom, str(pos)])

        if flag_with_id:
            value = line.split('\t')[2]
        else:
            value = 'T'

        if value not in db_dict[key]:
            db_dict[key].append(value)

    f.close()
    return db_dict


def load_genome(db):
    ''' load chromosome '''

    if db.endswith('.gz') or db.endswith('.gzip'):
        import gzip
        f = gzip.open(db, 'rb')
    else:
        f = open(db, 'r')

    db_dict = defaultdict(list)
    is_db_vcf = False

    for line in f:
        if line.startswith('#') or line == '':
            if 'fileformat=VCF' in line:
                is_db_vcf = True
            continue

        assert is_db_vcf == True, 'database must be in vcf format'

        chrom, pos = line.split('\t')[0:2]

        # resolve position when record is a deletion
        if len(line.split()[3]) > len(line.split()[4]):
            pos = resolve_db_position(args.input_type, int(pos))

        chrom.replace('chr', '')
        key = ','.join([chrom, str(pos)])

        if args.flag_with_id:
            value = line.split('\t')[2]
        else:
            value = 'T'

        if value not in db_dict[key]:
            db_dict[key].append(value)

    f.close()

    return db_dict


def load_db_positions(db_file, chrom):
    ''' load db positions into a dict '''

    if chrom is None:
        db_dict = load_genome(db_file)
    else:
        db_dict = load_chromosome(db_file, chrom)

    return db_dict


def write_header(infile, label, database, out, flag_with_id):
    ''' write the header for the output vcf '''

    db_hdr = '##{}_DB={}\n'.format(label, os.path.abspath(database))

    if flag_with_id:
        db_hdr += '##INFO=<ID={0},Number=.,Type=String,Description="{0} flag">\n'.format(label)
    else:
        db_hdr += '##INFO=<ID={0},Number=0,Type=Flag,Description="{0} flag">\n'.format(label)

    with open(infile, 'r') as f:

        # write the hdr and add new db descriptor
        hdr = []
        for line in f:
            if not line.startswith('#'):
                break
            hdr.append(line)

    if hdr:
        hdr = hdr[:-1] + [db_hdr] + [hdr[-1]]
    out.write(''.join(hdr))


# changed because "flagpos" seems confusing
# as some of the annotations are themselves
# 'flags'
def add_db_annotation(
        database, input_vcf, chromosome, output, label, flag_with_id, input_type
):
    ''' flag the vcf entries with
    information from the database annotator '''

    db_dict = None

    with open(output, 'w') as out:

        write_header(input_vcf, label, database, out, flag_with_id)

        with open(input_vcf, 'r') as f:
            for line in f:

                if line.startswith('#'):
                    continue

                line = line.rstrip().split('\t')
                chrom, pos = line[0:2]
                chrom = chrom.replace('chr', '')

                if chromosome and chrom != chromosome:
                    continue

                if not db_dict or not db_dict[0] == chrom:
                    db_dict = (chrom, load_chromosome(database, chromosome, input_type, flag_with_id))

                info = line[7]

                key = ','.join([chrom, pos])
                value = db_dict[1].get(key)

                if flag_with_id:
                    flag = ','.join(value) if value else '.'
                    info = '{};{}={}'.format(info, label, flag)
                else:
                    # according to vcf style, don't need anything
                    # on !value
                    if value:
                        info = '{};{}'.format(info, label)

                line[7] = info
                line = '\t'.join(line) + '\n'
                out.write(line)


def parse_args():
    chromosomes = map(str, range(1, 23))
    chromosomes.append('X')
    chromosomes.append('Y')

    parser = argparse.ArgumentParser(prog='flag positions',
                                     description='''flag a position in the INFO section of a vcf
                                                    file if it exists in a specified database''')

    parser.add_argument("--infile",
                        required=True,
                        help='''specify the path/to/VCF-formatted-input-file.vcf''')

    parser.add_argument("--input_type",
                        required=True,
                        choices=['indel', 'snv'],
                        help='''snv or indel''')

    parser.add_argument("--db",
                        required=True,
                        help='''specify the /path/to/db_file to be
                   used as a reference to determine whether
                   positions in input file is also in db''')

    parser.add_argument("--out",
                        required=True,
                        help='''specify the /path/to/out.vcf to save output to a file''')

    parser.add_argument("--label",
                        required=True,
                        help='''specify a label to be used e.g. "<label>=T"''')

    parser.add_argument("--flag_with_id",
                        action='store_true',
                        help='''use a list of ids as the flag for presence in the database
                           instead of 'T' ''')

    parser.add_argument("--chrom",
                        required=False,
                        choices=chromosomes,
                        default=None,
                        help='''chromosome''')

    args = parser.parse_args()

    args = vars(args)

    return args


if __name__ == '__main__':
    args = parse_args()
    add_db_annotation(
        args['db'], args['infile'], args['chrom'], args['out'],
        args['label'], args['flag_with_id'], args['input_type']
    )
