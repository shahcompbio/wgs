from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree


def build_interval_tree(data):
    '''
    load lumpy confidence intervals into tree
    '''
    itree = defaultdict(IntervalTree)

    for val in data:
        chrom, pos, ci, brk_id = val

        ci = ci.split(',')

        start = pos + float(ci[0])
        end = pos + float(ci[1]) + 1

        itree[chrom].addi(start, end, brk_id)

    return itree


def check_olp(interval_tree, chrom, pos):
    if interval_tree[chrom][pos]:
        return True
    else:
        return False


def load_data(infile):
    return pd.read_csv(infile)


def load_lumpy_into_tree(lumpy_df, confidence_interval=None):
    # Add lumpy breakpoint id to each zipped entry
    if confidence_interval:
        confidence_interval = '-{},{}'.format(confidence_interval, confidence_interval)
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, [confidence_interval] * len(lumpy_df)))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, [confidence_interval] * len(lumpy_df)))
    else:
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, lumpy_df['CIPOS']))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, lumpy_df['CIEND']))

    intervaltree = build_interval_tree(data)

    return intervaltree


def filter_destruct_on_lumpy(destruct, lumpy_tree):
    destruct['is_lumpy_filtered'] = False

    for idx in destruct.index:
        lumpy_brk_ids = dict()

        for break_end in ('1', '2'):
            chrom = destruct.loc[idx, 'chromosome_{break_end}'.format(break_end)]
            pos = destruct.loc[idx, 'position_{break_end}'.format(break_end)]

            for ival in lumpy_tree[chrom][pos]:
                lumpy_brk_ids[break_end] = set([ival[2] for ival in lumpy_tree[chrom][pos]])

        if len(lumpy_brk_ids['1'].intersection(lumpy_brk_ids['2'])) > 0:
            destruct.loc[idx, 'is_lumpy_filtered'] = True

    return destruct


def write(data, outfile):
    data.to_csv(outfile, index=False)


def consensus(destruct_infile, lumpy_infile, consensus, confidence_interval=None):
    destruct = load_data(destruct_infile)
    lumpy = load_data(lumpy_infile)

    lumpy = load_lumpy_into_tree(lumpy, confidence_interval=confidence_interval)

    destruct = filter_destruct_on_lumpy(destruct, lumpy)

    write(destruct, consensus)
