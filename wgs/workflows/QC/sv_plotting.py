import pandas as pd


def read(sv_data):
    svs = pd.read_csv(sv_data, sep=",")[['chromosome_1', 'position_1',
                                         'chromosome_2', 'position_2', 'rearrangement_type']]
    return add_colors_to_sv_data(svs)


def add_colors_to_sv_data(sv_data):
    types = ['foldback', 'unbalanced', 'duplication',
            'deletion', 'inversion', 'balanced']
    colors = [2, 3, 4, 1, 6, 8]
    sv_data["color"] = sv_data.rearrangement_type.replace(types, colors)
    return sv_data


def make_for_circos(sv_calls, outfile):
    svs = read(sv_calls)
    svs.to_csv(outfile, index=False, header=True, sep="\t")

##Example
#run make_for_circos() for filtered_consensus_calls sv output
#and use output in circos.R