import pandas as pd
import sys
import matplotlib as mp

def get_chroms(b, chrom = "chrom"):
    return b[chrom]

def get_positions(b, pos = "pos"):
    return b[pos]

def compare_csvs(breakpoints_1, breakpoints_2):
    b1 = pd.read_csv(breakpoints_1, sep = ",")
    b2 = pd.read_csv(breakpoints_2, sep = ",")

    b1 = b1[["chrom", "pos", "ref", "alt", "PR"]]
    b2 = b2[["chrom", "pos", "ref", "alt", "PR"]]
    
    for row1, row2 in zip(b1.iterrows(), b2.iterrows()):
        comparison = compare(row1, row2)
        print (row1[1]["PR"])
        print (comparison)

def main(breakpoints_1, breakpoints_2):
    b1 = pd.read_csv(breakpoints_1, sep = ",")
    b2 = pd.read_csv(breakpoints_2, sep = ",")

    chroms_b1 = get_chroms(b1)

if __name__ == "__main__":
    main()