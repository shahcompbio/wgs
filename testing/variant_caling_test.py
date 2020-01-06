import pandas as pd
import sys
import matplotlib.pyplot as plt
from wgs.workflows.variant_calling_consensus.scripts import vcfparser as parser
import vcf
from collections import defaultdict

def compare(call1, call2):
    '''
    compares two dictionaries key by key.
    Returns 0 if keys differ.
    :param call1: dictionary 1
    :param call2: dictionary 2
    '''
    score = 0
    diffs = {}

    # if row1 does not at least contain keys in row2
    if call1.keys() <= set(call2.keys()) :
        return -1

    for k,val in call2.items():
        val2 = call1[k]
        if val == val2: 
            score +=1
        else:
            diffs[k] = (val, val2)
            
    return score, diffs
    
def get_chroms_and_pos(variants, filt, reader):
    '''
    makes a dictionary of all the positions in a variant reader
    where the keys are chromosomes.
    :param variants: variant reader
    '''
    chroms_and_pos = defaultdict(set)
    
    for row in variants:
        if not row.startswith("#"):
            row = row.split()
            if len(row) >= 4 and filt(row[0],int(row[1]) - 1, int(row[1]), reader) == True: 
                chroms_and_pos[row[0]].add((int(row[1]) - 1, int(row[1])))
    return chroms_and_pos

def filt(chrom, start, end, reader):
    '''
    filter positions by location 
    and PR val
    :param chrom: chrom of call
    :param start: start of call
    :param end: end of call
    :param reader: call file reader obj
    '''
    if int(start) > 10000000:
        call = get_call(reader, chrom, start, end)
        c = next(call).INFO["PR"]
        if c > 0.55:
            return True
    else:
        return True

#just at 1 position ## assums 1 call per pos
def get_call(reader, chrom, start, end):
    '''
    get call at position chrom, start
    :param reader: reader obj
    :param chrom: chrom of call
    :param start: start of call
    :param end: end of call
    '''
    return reader.fetch(chrom = chrom, start = start, end = end)

def intersect_pos_sets(pos_set1, pos_set2, chroms):
    '''
    finds identical positions on chromosomes
    :param pos_set1: dict of chrom: positions
    :param pos_set2: dict of chrom: positions
    :param chroms: list of shared chromosomes
    '''
    return {chrom:pos_set1[chrom].intersection(pos_set2[chrom]) 
        for chrom in chroms} 
  
def difference_pos_sets(pos_set1, pos_set2, chroms):
    '''
    finds differing positions on chromosomes
    :param pos_set1: dict of chrom: positions
    :param pos_set2: dict of chrom: positions
    :param chroms: list of shared chromosomes
    '''

    return {chrom:pos_set1[chrom].symmetric_difference(pos_set2[chrom]) 
        for chrom in chroms}

def main(variants1, variants2):
    
    v1_parser = parser.VcfParser(variants1 + ".gz", None, None,None,None,None)
    v2_parser = parser.VcfParser(variants2 + ".gz",None,None,None,None,None)

    v1_reader= vcf.Reader(filename = variants1 + ".gz")
    v2_reader= vcf.Reader(filename = variants2 + ".gz")
    
    v1 = open(variants1)
    v2 = open(variants2)
    print ('Parsing data')
    pos_set1 = get_chroms_and_pos(v1, filt, v1_reader)
    pos_set2 = get_chroms_and_pos(v2, filt, v2_reader)
    print (pos_set1)
    print ('Checking call positions')

    if "Y" in pos_set1.keys():
        del pos_set1["Y"]
    if "Y" in pos_set2.keys():
        del pos_set2["Y"]

    if not pos_set1.keys() == pos_set2.keys():
        print ("different chroms {} {}".format(pos_set1.keys(), 
            pos_set2.keys()))
        raise Exception()
    
    chroms = pos_set1.keys()
    print (pos_set1- pos_set2)
    similarities = intersect_pos_sets(pos_set1, pos_set2, chroms)
    differences = difference_pos_sets(pos_set1, pos_set2, chroms)
    
    if differences:
        print ("different positions {}".format(differences))
        raise Exception()
    
    print ('Calls are at the same positions, checking call content')
    for chrom, positions in similarities.iteritems():
        for s,e in similarities[chrom]:
            record_v1 = v1_parser.parse_record(
                next(v1.fetch(chrom = chrom, start = s, end = e)))
            record_v2 = v2_parser.parse_record(
                next(v2.fetch(chrom = chrom, start = s, end = e)))
            compare(record_v1, record_v2)
    v1.close()
    v2.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
