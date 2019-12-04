import numpy as np
import pandas as pandas

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
