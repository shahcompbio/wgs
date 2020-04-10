from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree


def build_interval_tree(data):
    '''
    load lumpy confidence intervals into tree
    '''
    itree = defaultdict(IntervalTree)

    for val in data:
        chrom, pos, ci, brk_type = val

        ci = ci.split(',')

        start = pos + float(ci[0])
        end = pos + float(ci[1]) + 1

        itree[chrom].addi(start, end, brk_type)

    return itree


def check_olp(interval_tree, chrom, pos, destruct_strand):
    overlap = interval_tree[chrom][pos]

    overlapping_strands = [v[2] for v in overlap]

    if destruct_strand in overlapping_strands:
        return True

    return False



def load_data(infile):
    return pd.read_csv(infile, dtype={'chromosome_1': str, 'chromosome_2': str})


def load_lumpy_into_tree(lumpy_df, confidence_interval=None):
    # Add lumpy breakpoint id to each zipped entry
 
    if confidence_interval:
        confidence_interval = '-{},{}'.format(confidence_interval, confidence_interval)
        data = list(
            zip(lumpy_df.chromosome_1, lumpy_df.position_1, [confidence_interval] * len(lumpy_df), lumpy_df.strand_1))
        data += list(
            zip(lumpy_df.chromosome_2, lumpy_df.position_2, [confidence_interval] * len(lumpy_df), lumpy_df.strand_2))
    else:
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, lumpy_df['CIPOS'], lumpy_df['strand_1']))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, lumpy_df['CIEND'], lumpy_df['strand_2']))

    intervaltree = build_interval_tree(data)

    return intervaltree


def filter_destruct_on_lumpy(destruct, lumpy_tree):
    destruct = destruct[
        destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_1'], x['position_1'], x['strand_1']),
                       axis=1)]
    destruct = destruct[
        destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_2'], x['position_2'], x['strand_2']),
                       axis=1)]

    return destruct


def write(data, outfile):
    data.to_csv(outfile, index=False)


def consensus(destruct_infile, lumpy_infile, consensus, confidence_interval=None):
    destruct = load_data(destruct_infile)
    lumpy = load_data(lumpy_infile)

    lumpy = load_lumpy_into_tree(lumpy, confidence_interval=confidence_interval)

    destruct = filter_destruct_on_lumpy(destruct, lumpy)

    write(destruct, consensus)


