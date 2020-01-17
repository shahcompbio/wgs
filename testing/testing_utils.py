import pandas as pd
import sys
import os.path
import yaml
import numpy as np

def read_in_files(files, yamls, seps):
    '''
    gets assumed files from directory dir
    :param dir: path to titan files
    '''
    if isinstance(files, list):
        return [read_in(f, yaml, sep ) for yaml, f, sep in zip(yamls, files, seps)]
    if isinstance(files, dict):
        return {name: read_in(f, yamls[name], seps[name]) for name, f in files.items()}

def get_index(df1, df2, on):
    return df1[on].stack().isin(df2[on].stack().values).unstack()

def df_difference(df1, df2, on):
    index1 = ~get_index(df1, df2, on)
    index2 = ~get_index(df2, df1, on)

    return df1[index1[on[0]] | index1[on[1]]].dropna(), \
        df2[index2[on[0]] | index2[on[1]]].dropna()

def df_intersection(df1, df2, on):
    index1 = get_index(df1, df2, on)
    index2 = get_index(df2, df1, on)

    return df1[index1[on[0]] & index1[on[1]]].dropna(), \
        df2[index2[on[0]] & index2[on[1]]].dropna()

def read_in(file, yaml_dtypes=None, sep = "," ):
    '''
    parse in csv with yaml dtypes
    :param file: csv with csv.yaml 
    '''
    if yaml_dtypes:
        dtypes = {col['name']:col['dtype'] 
            for col in yaml.load(open(yaml_dtypes))["columns"]}
        return pd.read_csv(file, dtype = dtypes, sep = sep)
    else:
        return pd.read_csv(file, sep = sep)

def check_cols(df1, df2):
    return df1.columns.tolist() == df2.columns.tolist()


def compare_sames(df1, df2):
    '''
    get col-wise difference between two
    pandas dataframes
    :param df1: pandas dataframe
    :param df2: pandas dataframe
    '''
    if df1 == df2:
        return True
    differences = {}
    if df1.columns == df2.columns:
        differences = {col:{} for col in df1.columns}
        for (name1, col1), (name2, col2) in zip(df1.items(), df2.items()): 
            if col1 != col2:
                 for (i,row1), (i2,row2) in zip(enumerate(col1), enumerate(col2)):
                     if row1 != row2:
                         differences[col1].add((i, row1, row2))      
    return False, differences