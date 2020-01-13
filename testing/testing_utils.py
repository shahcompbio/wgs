import pandas as pd
import sys
import os.path
import yaml

def read_in_files(files, yamls):
    '''
    gets assumed files from directory dir
    :param dir: path to titan files
    '''
    if isinstance(files, list):
        return [read_in(f, yaml ) for yaml, f in zip(yamls, files)]
    if isinstance(files, dict):
        return {name: read_in(f, yamls[name]) for name, f in files.items()}


def read_in(file, yaml_dtypes=None):
    '''
    parse in csv with yaml dtypes
    :param file: csv with csv.yaml
    '''
    if yaml_dtypes:
        dtypes = {col['name']:col['dtype'] 
            for col in yaml.load(open(yaml_dtypes))["columns"]}
        return pd.read_csv(file, dtype = dtypes, sep = "\t")
    else:
        return pd.read_csv(file, sep = "\t")


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
        for (name1, col1), (name2, col2) in zip(df1.iteritems(), df2.iteritems()): 
            if col1 != col2:
                 for (i,row1), (i2,row2) in zip(enumerate(col1), enumerate(col2)):
                     if row1 != row2:
                         differences[col1].add((i, row1, row2))      
    return False, differences