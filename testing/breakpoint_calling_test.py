import pandas as pd
import sys
import matplotlib as mp
import testing_utils as utils 
import os.path


def df_symmetric_difference(df1, df2, columns):
    '''
    gets the difference between dataframe 1 and 2 on given columns
    :param df1: pandas df
    :param df2: pandas df
    :columns: shared columns to compare between df1, df2
    '''
    return pd.concat([df1, df2]).drop_duplicates(keep = False, subset = columns)

def df_sames(df1, df2, columns):
    '''
    gets the same records between df1 and df2
    :param df1: pandas df
    :param df2: pandas df
    :columns: shared columns to compare between df1, df2
    '''
    both= pd.concat([df1, df2])
    return both[both.duplicated(keep = False, subset = columns)]


def main(dir1, dir2, files1, files2):
    files1 = {name: os.path.join(dir1, f) for name, f in files1.items()}
    files2 = {name: os.path.join(dir2, f) for name, f in files2.items()}
 
    yamls = {name:f+ ".csv.gz.yaml" for name, f in files1.items()}

    dfs1 = utils.read_in_files(files1, yamls)
    dfs2 = utils.read_in_files(files2, yamls)

    breakpoints1 = dfs1["breakpoints"]
    breakpoints2 = dfs2["breakpoints"]

    lib1 = dfs1["library"]
    lib2 = dfs2["library"]

    assert breakpoints1.columns.tolist() == breakpoints2.columns.tolist()
    breakpoint_cols = breakpoints1.columns

    assert lib1.columns.tolist() == lib2.columns.tolist()
    lib_cols = lib1.columns

    breakpoints1 = breakpoints1.sort_values(
        ["chromosome_1", "chromosome_2","position_1", "position_2"], 
        ascending =[True, True, True, True]
    )
    # breakpoints1.to_csv("bp1_sorted.tsv",index = False, header = True, sep = "\t")
    breakpoints2 = breakpoints2.sort_values(
        ["chromosome_1","chromosome_2", "position_1", "position_2"], 
        ascending =[True, True, True, True]
    )
    # breakpoints2.to_csv("bp2_sorted.tsv",index = False, header = True, sep = "\t")
    print ("diff")

    pos_set1 = {tuple(x) for x in breakpoints1[["chromosome_1","chromosome_2", "position_1", "position_2"]].values}
    pos_set2 = {tuple(x) for x in breakpoint:qs2[["chromosome_1","chromosome_2", "position_1", "position_2"]].values}

    d1 = pos_set1-pos_set2
    d2 = pos_set2-pos_set1

    d1 = pd.DataFrame(list(d1), columns=["chromosome_1","chromosome_2", "position_1", "position_2"])
    d2 = pd.DataFrame(list(d2), columns=["chromosome_1","chromosome_2", "position_1", "position_2"])

    d1.to_csv("posdiffs1.tsv",index = False, header = True, sep = "\t")
    d2.to_csv("posdiffs2.tsv",index = False, header = True, sep = "\t")

    sames = df_sames(breakpoints1, breakpoints2, columns = ["chromosome_1","chromosome_2", "position_1", "position_2"])
    sames.to_csv("sames.tsv",index = False, header = True, sep = "\t")
    # diffs = df_symmetric_difference(breakpoints1[["chromosome_1","chromosome_2", "position_1", "position_2"]], breakpoints2[["chromosome_1","chromosome_2", "position_1", "position_2"]])
    # diffs.to_csv("posdiffs.tsv",index = False, header = True, sep = "\t")
    
    return
    # for col in lib_cols:
    #     col1 = lib1[col]
    #     col2 = lib2[col]
    #     assert col1 == col2




if __name__ == "__main__":
    main(dir1 = "/juno/work/shah/dgrewal/work/COMP_SCP/bkp/temp/sample_id/SA123/destruct",
        dir2 = "/juno/work/shah/dgrewal/work/COMP_SCP/kronos/kronos_results/SA1229T_destruct",
        files1 = {"breakpoints": "raw_breakpoints", "library": "raw_library"},
        files2 = {"breakpoints": "TASK_1_RUN_breakpoint_table.txt", "library": "TASK_1_RUN_breakpoint_library_table.txt"})