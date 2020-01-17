import pandas as pd
import sys
import os.path
import yaml
from functools import reduce
import testing_utils as utils
import numpy as np
import matplotlib

matplotlib.use('agg')

import matplotlib.pyplot as plt
def get_dfs(dir):
    '''
    gets assumed files from directory dir
    :param dir: path to titan files
    '''
    files = dir["files"]
    path  = dir["path"]
    seps  = {} 
    yamls = {}

    for label, f in files.items():
        files[label] = os.path.join(path, f[0])
        seps[label]  = f[1]
        yamls[label] = files[label] + ".yaml"

    return utils.read_in_files(files, yamls, seps)
     
def plot_differences(df1, df2):
    '''
    plot a histogram for the differences of every column
    between df1 and df2. should have same chrom / pos
    :param df1: pandas dataframe
    :param df2: pandas dataframe
    '''
    #reset indexes
    df1.reset_index(inplace = True) 
    df2.reset_index(inplace = True) 

    #only take numerics
    df1 = df1.select_dtypes(include=np.number)
    df2 = df2.select_dtypes(include=np.number)

    diffs = df1.subtract(df2).abs()
    means = (df1.mean(axis=0) + df2.mean(axis=0))/2
    diffs = diffs.div(means, axis = 1)

    fig, axes = plt.subplots(len(diffs.columns), 1, figsize=(10, 40))    
    diffs.hist( bins = 20, ax=axes)

    return fig 

def plot_ref_count_differers(df1, df2):
    '''
    plot only the differing entries for the ref count 
    fields.
    :param df1: pandas dataframe
    :aram df2: pandas dataframe
    '''
    df1_rc_diff, df2_rc_diff = utils.df_difference(df1, df2, 
        ["NRefCount", "RefCount"])
    return plot_differences(df1_rc_diff, df2_rc_diff)


def merge_on_titan(dfs):
    '''
    merge dataframs in dfs on 
    TITAN_state and TITAN_call.
    :param dfs: list of pandas datafram objs
    '''
    return reduce(lambda left,right: 
        pd.merge(left,right,on=["TITAN_state","TITAN_call"]), dfs)

def main(kronos, wgs_new):
    kronos_dfs = get_dfs(kronos)
    wgs_new_dfs = get_dfs(wgs_new)

    kronos_segs = kronos_dfs["segs"]
    wgs_new_segs = wgs_new_dfs["segs"]

    kronos_markers = kronos_dfs["markers"]
    wgs_new_markers = wgs_new_dfs["markers"]

    #these are the overlapping columns 
    kronos_parsed = kronos_dfs["parsed"][["chromosome", "start", "stop"]]
    wgs_new_parsed = wgs_new_dfs["parsed"][["chrom", "start", "end"]]

    if utils.check_cols(kronos_segs,wgs_new_segs):
        segs_kronos_diff, segs_wgs_new_diff = utils.df_difference(
            kronos_segs, wgs_new_segs, 
            on=['Chromosome', 'Start_Position(bp)']
        )
        segs_kronos_same, segs_wgs_new_same = utils.df_intersection(
            kronos_segs, wgs_new_segs, 
            on=['Chromosome', 'Start_Position(bp)']
        )

        fig = plot_differences(segs_kronos_same, segs_wgs_new_same)
        fig.savefig("segments_field-wise_differences.pdf")
        
        segs_kronos_diff.to_csv("differing_segments_kronos.csv", sep = ",")
        segs_wgs_new_diff.to_csv("differing_segments_wgs_new.csv", sep = ",")

    if utils.check_cols(kronos_markers, wgs_new_markers):
        markers_kronos_diff, markers_wgs_new_diff = utils.df_difference(
            kronos_markers, wgs_new_markers, 
            on = ['Chr', 'Position']
        )
        markers_kronos_same, markers_wgs_new_same = utils.df_intersection(
            kronos_markers, wgs_new_markers, 
            on = ['Chr', 'Position']
        )

        fig = plot_differences(markers_kronos_same, markers_wgs_new_same)
        fig.savefig("markers_field-wise_differences.pdf")
        
        fig = plot_ref_count_differers(markers_kronos_same, markers_wgs_new_same)
        fig.savefig("markers_differences_readcount.pdf")

        markers_kronos_diff.to_csv("differing_markers_kronos.csv", sep = ",")
        markers_wgs_new_diff.to_csv("differing_markers_wgs_new.csv", sep = ",")

    if utils.check_cols(kronos_parsed,wgs_new_parsed):
        parsed_kronos_diff, parsed_wgs_new_diff = utils.df_difference(
            kronos_parsed, wgs_new_parsed, 
            on=['Chromosome', 'Start_Position(bp)']
        )
        parsed_kronos_same, parsed_wgs_new_same = utils.df_intersection(
            kronos_parsed, wgs_new_parsed,
            on=['Chromosome', 'Start_Position(bp)']
        )

        fig = plot_differences(parsed_kronos_same, parsed_wgs_new_same)
        fig.savefig("segments_field-wise_differences.pdf")

        parsed_kronos_diff.to_csv("differing_parsed_kronos.csv", sep = ",")
        parsed_wgs_new_diff.to_csv("differing_parsed_wgs_new.csv", sep = ",")
 
if __name__ == "__main__":
    KRONOS = {
        "path": "/juno/work/shah/abramsd/SA500titan/output/copynumber/SA500/titan/optimal",
        "files":{
        "segs": ["kronos_segs.tsv", "\t"], 
        "markers": ["kronos_markers.tsv", "\t"], 
        "parsed": ["TASK_12__5_titan_parsed.tsv", "\t"]}}
    WGS_NEW = {
        "path": "/juno/work/shah/abramsd/SA500titan/output/copynumber/SA500/titan",
        "files":{
        "segs":["SA500_titan_segs.csv.gz", "\t"], 
        "markers": ["SA500_titan_markers.csv.gz", "\t"] , 
        "parsed":["SA500_titan_parsed.csv.gz", ","]}}
    main(KRONOS, WGS_NEW)