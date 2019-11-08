import numpy as np
import pandas as pd


def load_destruct_calls(infile):
    return pd.read_csv(infile, sep='\t')


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


def parser(infile, outfile, foldback_threshold=30000):
    df = load_destruct_calls(infile)

    df = fix_dgv_ids(df)

    df = add_break_distance(df)

    df = reclassify_inversions(df)

    df = classify_foldbacks(df, foldback_threshold)

    write(df, outfile)
