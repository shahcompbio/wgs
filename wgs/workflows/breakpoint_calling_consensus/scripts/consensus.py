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


def compare_lumpy_destruct_type(lumpy_type, destruct_type):
    if lumpy_type not in ['DEL', 'INV', 'DUP', 'BND']:
        raise Exception('Unknown lumpy type: {}'.format(lumpy_type))

    if destruct_type not in ['deletion', 'inversion', 'duplication', 'translocation']:
        raise Exception('Unknown destruct type: {}'.format(destruct_type))

    if lumpy_type == 'DEL' and destruct_type == 'deletion':
        return True
    elif lumpy_type == 'INV' and destruct_type == 'inversion':
        return True
    elif lumpy_type == 'DUP' and destruct_type == 'duplication':
        return True
    elif lumpy_type == 'BND' and destruct_type == 'translocation':
        return True
    else:
        return False


def check_olp(interval_tree, chrom, pos, destruct_type):
    overlap = interval_tree[chrom][pos]

    lumpy_types = [v[2] for v in overlap]
    lumpy_types = [compare_lumpy_destruct_type(v, destruct_type) for v in lumpy_types]

    if lumpy_types and any(lumpy_types):
        return True
    else:
        return False


def load_data(infile):
    return pd.read_csv(infile)


def load_lumpy_into_tree(lumpy_df, confidence_interval=None):
    # Add lumpy breakpoint id to each zipped entry

    if confidence_interval:
        confidence_interval = '-{},{}'.format(confidence_interval, confidence_interval)
        data = list(
            zip(lumpy_df.chromosome_1, lumpy_df.position_1, [confidence_interval] * len(lumpy_df), lumpy_df.SVTYPE))
        data += list(
            zip(lumpy_df.chromosome_2, lumpy_df.position_2, [confidence_interval] * len(lumpy_df), lumpy_df.SVTYPE))
    else:
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, lumpy_df['CIPOS'], lumpy_df['SVTYPE']))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, lumpy_df['CIEND'], lumpy_df['SVTYPE']))

    intervaltree = build_interval_tree(data)

    return intervaltree


def filter_destruct_on_lumpy(destruct, lumpy_tree):
    destruct = destruct[
        destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_1'], x['position_1'], x['rearrangement_type']),
                       axis=1)]
    destruct = destruct[
        destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_2'], x['position_2'], x['rearrangement_type']),
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
