import numpy as np
import pandas as pd


def load_blacklist(blacklist):
    blacklist = pd.read_csv(blacklist, sep='\t')
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


def annotate_low_mappability(destruct_calls, blacklist):
    '''
    adds low mappability annotation to destruct
    calls.
    :param destruct_calls: pandas df of destruct
    calls
    '''
    is_low_mapp_1 = is_low_mappability(destruct_calls,
                                       blacklist, "chromosome_1", "position_1")

    is_low_mapp_2 = is_low_mappability(destruct_calls,
                                       blacklist, "chromosome_2", "position_2")

    low_mapp_indexes = list(set(is_low_mapp_1) | set(is_low_mapp_2))

    low_mappability_annotation = generate_low_mappability_annotation(
        low_mapp_indexes, len(destruct_calls.index))

    destruct_calls["is_low_mappability"] = low_mappability_annotation

    return destruct_calls


def main(infile, outfile, mappability_blacklist):
    blacklist = load_blacklist(mappability_blacklist)

    brk = pd.read_csv(
        infile, sep='\t',
        converters={'chromosome_1': str, 'chromosome_2': str}
    )

    brk = annotate_low_mappability(brk, blacklist)

    brk.to_csv(outfile, sep='\t', index=False, header=True)
