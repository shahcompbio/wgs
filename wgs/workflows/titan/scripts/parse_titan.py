from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree


def build_interval_tree(data):
    '''
    load lumpy confidence intervals into tree
    '''
    itree = defaultdict(IntervalTree)

    data = zip(data['chrom'], data['start'], data['end'])

    for val in data:
        chrom, start, end = val

        if start == end:
            continue

        itree[chrom].addi(start, end)

    return itree


def load_segment_data(infile):
    df = pd.read_csv(infile, sep='\t', dtype={'Chromosome': str})

    df['chrom'] = df['Chromosome']
    df['start'] = df['Start_Position(bp)']
    df['end'] = df['End_Position(bp)']
    df['width'] = df['Length(bp)']

    del df['Pygenes(gene_id,gene_name;)']
    del df['Chromosome']
    del df['Start_Position(bp)']
    del df['End_Position(bp)']

    return df


def load_markers_data(infile):
    return pd.read_csv(infile, sep='\t', dtype={'Chr': str})


def get_marker_counts(markers, segs):
    segs_tree = build_interval_tree(segs)

    counts_dict = {}

    for (chrom, pos) in zip(markers.Chr, markers.Position):

        ival = segs_tree[chrom][pos]

        if not list(ival):
            continue

        start, stop, _ = list(ival)[0]

        if (chrom, start, stop) not in counts_dict:
            counts_dict[(chrom, start, stop)] = 0

        counts_dict[(chrom, start, stop)] += 1

    return counts_dict


def add_counts_to_segs(segs, marker_counts):
    segs["count_markers"] = segs.apply(lambda x: marker_counts.get((x['chrom'], x['start'], x['end']), 0), axis=1)
    return segs


def write(df, output):
    df.to_csv(output, index=False)


def parser(markers, segs, output):
    segs = load_segment_data(segs)
    markers = load_markers_data(markers)

    marker_counts = get_marker_counts(markers, segs)

    add_counts_to_segs(segs, marker_counts)

    write(segs, output)


if __name__ == "__main__":
    parser('titan_markers.txt', 'titan_segs.txt', 'titan_parsed.txt')
