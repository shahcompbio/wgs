import numpy as np
import pandas as pd
from wgs.utils import helpers


def load_blacklist(blacklist):
    blacklist = pd.read_csv(blacklist, sep='\t', dtype={'chromosome': str})
    return blacklist


def in_any_region(pos, blacklist):
    '''
    checks if pos are in any of the regions
    covered by blacklist. Expects that all positions
    in both fall on same chromosome.
    :param pos: pandas dataframe of blacklisted regions
    :param pos: integer of pos to be searched in blacklist
    '''
    if not len(blacklist):
        return []

    idx = np.searchsorted(blacklist['start'], pos, side='right') - 1
    return pos <= blacklist.loc[blacklist.index[idx], 'end'].values


def is_low_mappability(muts, blacklist, chromname, posname):
    '''
    gets entries in muts that fall within low mappability regions
    of blacklist.

    :param muts: csv with variant calls
    :param blacklist: csv of blacklist locations
    :param chromname: column name for chromosomes in muts
    :param posname: column name for positions in muts
    '''
    low_mappability_calls = []
    for chrom in muts[chromname].unique():
        chrom_muts = muts[muts[chromname] == chrom]
        chrom_blacklist = blacklist[blacklist['chromosome'] == str(chrom)]
        is_low_mapp = in_any_region(chrom_muts[posname], chrom_blacklist)
        low_mappability_calls.extend(chrom_muts.index[is_low_mapp])
    return low_mappability_calls


def generate_low_mappability_annotation(low_mappability_indexes, size):
    '''
    generates a list of booleans
    where the entries at low_mappability_indexes
    are True

    :param low_mappability_indexes: indexes to make True, in low mappability
    regions
    :param size: size of list
    '''
    return [True if i in low_mappability_indexes else False for i in range(size)]


def load_vcf_file(vcf_file):
    vcf_data = []
    with helpers.GetFileHandle(vcf_file) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue

            line = line.strip().split()

            vcf_data.append(line)

            if len(vcf_data) > 1e6:
                yield vcf_data
                vcf_data = []
    yield vcf_data


def get_vcf_header(vcf_file):
    vcf_data = []
    with helpers.GetFileHandle(vcf_file) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                vcf_data.append(line)
                continue
            break
    return vcf_data


def update_vcf_header(header, db_file):
    assert header[-1][0] == '#' and header[-1][1] != '#'

    updated_header = header[:-1]

    updated_header.append('##LOW_MAPPABILITY_DB={}'.format(db_file))
    updated_header.append('##INFO=<ID=LOW_MAPPABILITY,Number=0,Type=Flag,Description="low mappability position">')

    updated_header.append(header[-1])

    return updated_header


def write_to_file(output, lines):
    for line in lines:
        assert isinstance(line, str)
        line = line.strip('\n')
        output.write(line + '\n')


def annotate_vcf_data(vcf_data, blacklist):
    chroms = []
    positions = []
    for line in vcf_data:
        chroms.append(line[0])
        positions.append(int(line[1]))

    call_locations = pd.DataFrame({"chromosome": chroms, "positions": positions})

    low_mapp_indexes = is_low_mappability(call_locations, blacklist, "chromosome", "positions")

    annotated_data = []
    for i, line in enumerate(vcf_data):

        annotation = True if i in low_mapp_indexes else False
        if annotation:
            line[7] += ';LOW_MAPPABILITY'

        line = '\t'.join(line)
        annotated_data.append(line)

    return annotated_data


def main(infile, outfile, mappability_blacklist):
    blacklist = load_blacklist(mappability_blacklist)

    vcf_header = get_vcf_header(infile)
    vcf_header = update_vcf_header(vcf_header, mappability_blacklist)

    with helpers.GetFileHandle(outfile, 'wt') as vcf_writer:
        write_to_file(vcf_writer, vcf_header)

        for vcf_data in load_vcf_file(infile):
            annotated_vcf_data = annotate_vcf_data(vcf_data, blacklist)

            write_to_file(vcf_writer, annotated_vcf_data)
