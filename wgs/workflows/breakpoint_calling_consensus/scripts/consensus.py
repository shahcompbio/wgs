from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree


def build_interval_tree(data):
    '''
    load lumpy confidence intervals into tree
    '''
    itree = defaultdict(IntervalTree)

    for val in data:
        chrom, pos, ci = val

        ci = ci.split(',')

        start = pos + float(ci[0])
        end = pos + float(ci[1]) + 1

        itree[chrom].addi(start, end)

    return itree


def check_olp(interval_tree, chrom, pos):
    if interval_tree[chrom][pos]:
        return True
    else:
        return False


def load_data(infile):
    return pd.read_csv(infile)


def load_lumpy_into_tree(lumpy_df, confidence_interval=None):
    if confidence_interval:
        confidence_interval = '-{},{}'.format(confidence_interval, confidence_interval)
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, [confidence_interval]*len(lumpy_df)))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, [confidence_interval]*len(lumpy_df)))
    else:
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, lumpy_df['CIPOS']))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, lumpy_df['CIEND']))

    intervaltree = build_interval_tree(data)

    return intervaltree


def filter_destruct_on_lumpy(destruct, lumpy_tree):
    destruct = destruct[destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_1'], x['position_1']), axis=1)]
    destruct = destruct[destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_2'], x['position_2']), axis=1)]

    return destruct


def write(data, outfile):
    data.to_csv(outfile, index=False)


def consensus(destruct_infile, lumpy_infile, consensus, confidence_interval=None):
    destruct = load_data(destruct_infile)
    lumpy = load_data(lumpy_infile)

    lumpy = load_lumpy_into_tree(lumpy, confidence_interval=confidence_interval)

    destruct = filter_destruct_on_lumpy(destruct, lumpy)

    write(destruct, consensus)
