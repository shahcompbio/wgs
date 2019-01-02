import pandas as pd

def filter_destruct_breakpoints(breakpoints, filtered_breakpoints, min_num_reads):
    breakpoints = pd.read_csv(breakpoints, sep='\t')

    breakpoints = breakpoints[breakpoints['num_reads'] >= min_num_reads]

    breakpoints = breakpoints[[
        'prediction_id',
        'chromosome_1',
        'strand_1',
        'position_1',
        'chromosome_2',
        'strand_2',
        'position_2',
    ]]

    breakpoints.to_csv(filtered_breakpoints, sep='\t', index=False)
