'''
@author:  Celia Siu
created: 2014-07-31 (Thursday)
last_updated: 12 Feb 16 by dgrewal

'''
import argparse
import os
from collections import defaultdict
import warnings

version = '1.3.1'


def resolve_db_position(input_type, pos):
    '''
    adjust databse position to match convention of
    input file type for indexing a deletion
    '''
    if input_type == 'snv':
        pos += 1

    return pos


def load_chromosome(db, chromosome):
    ''' load genome reference '''
    import pysam
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
            pos = resolve_db_position(args.input_type, int(pos))

        chrom.replace('chr','')
        key = ','.join([chrom, str(pos)])

        if args.flag_with_id:
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

    db_dict   = defaultdict(list)
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

        chrom.replace('chr','')
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


def write_header(out):
    ''' write the header for the output vcf '''

    db_hdr =  '##{}_DB={}\n'.format(args.label, os.path.abspath(args.db))
    db_hdr += '##INFO=<ID={0},Number=.,Type=String,Description="{0} flag">\n'.format(args.label)

    with open(args.infile, 'r') as f:

        # write the hdr and add new db descriptor
        hdr = []
        for line in f:
            if not line.startswith('#'):
                break
            hdr.append(line)

    if hdr:
        hdr = hdr[:-1] + [db_hdr] + [hdr[-1]]
    out.write(''.join(hdr))


def flag_positions(chromosome, db_file):
    ''' flag the positions '''

    db_dict  =  None

    with open(args.out, 'w') as out:

        write_header(out)

        with open(args.infile, 'r') as f:
            for line in f:

                if line.startswith('#'):
                    continue

                line = line.rstrip().split('\t')
                chrom, pos = line[0:2]
                chrom = chrom.replace('chr','')

                if chromosome and chrom != chromosome:
                    continue

                if not db_dict or not db_dict[0] == chrom:
                    db_dict = (chrom, load_chromosome(db_file, chrom))

                info = line[7]

                key   = ','.join([chrom, pos])
                value = db_dict[1].get(key)

                if value:
                    flag = '[' + ','.join(value) + ']'
                else:
                    flag = 'F'

                info = '{};{}={}'.format(info, args.label, flag)
                line[7] = info
                line = '\t'.join(line)+'\n'
                out.write(line)


def main():
    ''' main function '''

#    db_dict = load_db_positions(args.db, args.chrom)
    flag_positions(args.chrom, args.db)


if __name__ == '__main__':

    chromosomes = map(str,range(1,23))
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
                   choices=['indel','snv'],
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


    args, unknown = parser.parse_known_args()

    main()
