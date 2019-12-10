import numpy as np
import pandas as pd


def load_destruct_calls(infile):
    return pd.read_csv(infile, sep='\t', dtype={'chromosome_1': str, 'chromosome_2': str})


def fix_dgv_ids(df):
    df['dgv_ids'] = df['dgv_ids'].replace(' ', '')
    return df


def add_break_distance(df, ):
    df['break_distance'] = np.absolute(df['position_1'] - df['position_2'])

    df['break_distance'][df['chromosome_1'] != df['chromosome_2']] = np.inf

    return df


def reclassify_inversions(df):
    df.loc[(df['type'] == 'inversion') & (df['rearrangement_type'] == 'unbalanced'), 'rearrangement_type'] = 'inversion'
    return df


def classify_foldbacks(df, threshold):
    df['intra'] = 0

    df.loc[(df['chromosome_1'] == df['chromosome_2']), 'intra'] = 1

    df.loc[
        (df['type'] == 'inversion') &
        (df['rearrangement_type'] == 'inversion') &
        (df['break_distance'] <= threshold) &
        (df['intra'] == 1), 'rearrangement_type'] = 'foldback'
    df.loc[
        (df['type'] == 'inversion') &
        (df['rearrangement_type'] == 'unbalanced') &
        (df['break_distance'] <= threshold) &
        (df['intra'] == 1), 'rearrangement_type'] = 'foldback'

    df = df.drop('intra', 1)

    return df


def write(df, outfile):
    df.to_csv(outfile, index=False)


def filter_calls(data, filters):
    if filters['readsupport_threshold']:
        data = data[data['num_reads'] > filters['readsupport_threshold']]

    if filters['chromosomes']:
        data = data[data['chromosome_1'].isin(filters['chromosomes'])]
        data = data[data['chromosome_2'].isin(filters['chromosomes'])]

    if filters['deletion_size_threshold']:
        data = data[
            (data['rearrangement_type'] != 'deletion') | (data['break_distance'] >= filters['deletion_size_threshold'])]

        data = data[
            (data['type'] != 'deletion') | (data['break_distance'] >= filters['deletion_size_threshold'])]

    if filters['break_distance_threshold']:
        data = data[data['break_distance'] > filters['break_distance_threshold']]

    return data


def parser(infile, outfile, filters, foldback_threshold=30000):
    df = load_destruct_calls(infile)

    df = fix_dgv_ids(df)

    df = add_break_distance(df)

    df = reclassify_inversions(df)

    df = classify_foldbacks(df, foldback_threshold)

    df = filter_calls(df, filters)

    write(df, outfile)
