import random
import vcf
import gzip
import os
import consensus
import pypeliner
from shutil import copyfile
import csv
from collections import namedtuple
import pandas as pd

def _get_test_record():
    '''
    make a dummy, correctly formatted snv with correct data
    Parameters
    ---------
    Returns
    -------
    '''
    test_record = [1, '.', "a", 'b', 34, 'a',1,2,3]
    chrom = random.choice(list(map(str, list(range(1, 23)) + ["X"])))
    pos = random.randint(0, 100000)
    alt = random.choice(["A", "T", "C", "G"])
    ref = random.choice(["A", "T", "C", "G"])
    return chrom, pos, ref, alt, test_record


def _get_indel():
    '''
    make a dummy, correctly formatted snv with correct data
    Parameters
    ---------
    Returns
    -------
    '''

    chrom = random.choice(list(map(str, list(range(1, 23)) + ["X"])))
    pos = random.randint(0, 100000)
    id_count = random.randint(0, 100000)
    test_record = [[ ".", ".", 1, 2, 3, 4, 5, 6, id_count], "T", "TC"]
    return chrom, pos, id_count, test_record


def _get_test_model_call():
    #make example  call
    Call = namedtuple("Call", "DP RO AO RC AC AD TAR TIR TU NoneU")

    normal_data = vcf.model._Call(None, "normal", Call(1,1,2,1,2, [1,2], [3],[4], [5], [6]) )
    tumor_data = vcf.model._Call(None, "tumor", Call(1,1,2,1,2, [1,2], [3],[4],[5], [6]) )

    freebayes_museq_rtg_example = vcf.model._Record("1", 10, 1, "T", [None], 1, 1, 
        None, None, None, [normal_data, tumor_data])

    samtools_example = vcf.model._Record("1", 10, 1, "T", [None], 1, 1, {"DP":1}, 
        None, None, [normal_data, tumor_data])

    return freebayes_museq_rtg_example, samtools_example

##################################
# test consensus.snv_consensus() #
##################################

def test_snv_consensus_case_1():
    '''
    test  consensus.snv_consensus() with identical record between museq/strelka
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    strelka = test_record
    mutect =[]

    consensus_data = consensus.snv_consensus(museq, strelka,mutect)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_count", "qual", "filter", 
        "tr", "ta", "td", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_snv_consensus_case_2():
    '''
    test  consensus.snv_consensus() with identical record between museq/mutect
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    strelka = []
    mutect = test_record

    consensus_data = consensus.snv_consensus(museq, strelka,mutect)
    
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_count", "qual", "filter", 
        "tr", "ta", "td", "nr", "na", "nd"])    

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_snv_consensus_case_3():
    '''
    test  consensus.snv_consensus() with identical record between strelka/mutect
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = []
    strelka = test_record
    mutect = test_record

    consensus_data = consensus.snv_consensus(museq, strelka,mutect)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_count", "qual", "filter", 
        "tr", "ta", "td", "nr", "na", "nd"])    

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_snv_consensus_case_4():
    '''
    test  consensus.snv_consensus() with no identical records
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    museq = {(chrom, pos + 1, ref, alt): record}

    strelka = {(chrom, pos + 2, ref, alt): record}

    mutect = {(chrom, pos + 3, ref, alt): record}


    consensus_data = consensus.snv_consensus(museq, strelka, mutect)

    assert consensus_data == []


def test_snv_consensus_case_5():
    '''
    test  consensus.snv_consensus() with empty data
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_data = consensus.snv_consensus([], [], [])

    assert consensus_data == []


####################################
# test consensus.indel_consensus() #
####################################


def test_indel_consensus_case_1():
    '''
    test consensus.indel_consensus() with identical record between strelka/mutect
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, id_count, record = _get_indel()

    test_record = {(chrom, pos, id_count): record}

    strelka = test_record
    mutect = test_record

    consensus_data = consensus.indel_consensus(strelka, mutect)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_count", "qual", "filter", 
        "tr", "ta", "td", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom":"str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_indel_consensus_case_2():
    '''
    test  consensus.snv_consensus() with empty data
    Parameters
    testdir: str directory
    ---------
    Returns
    -------
    '''
    consensus_data = consensus.indel_consensus([], [])

    assert consensus_data == []


def test_indel_consensus_case_3():
    '''
    test  consensus.snv_consensus() with no shared records, should output both
    Parameters
    testdir: str directory
    ---------
    Returns
    -------
    '''
    chrom, pos, id_count, record = _get_indel()

    strelka = {(chrom, pos + 2, id_count): record}

    mutect = {(chrom, pos + 3, id_count): record}

    consensus_data = consensus.indel_consensus(strelka, mutect)

    assert len(consensus_data) == 2


def test_indel_consensus_case_4():
    '''
    test consensus.indel_consensus() with no shared records
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, id_count, record = _get_indel()

    test_record = {(chrom, pos, id_count): record}

    strelka = []
    mutect = test_record

    consensus_data = consensus.indel_consensus(strelka, mutect)

    assert len(consensus_data) == 1


####################################
# test consensus.normalize() #
####################################

def test_normalization_case_1():
    '''
    test consensus.normalize() with identical ref/alt 
    expectation: no changes in  normalization
    Parameters
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize("T", "T")
    assert consensus_ref == "T"
    assert consensus_alt == "T"


def test_normalization_case_2():
    '''
    test consensus.normalize() with differing ref/alt length 1
    expectation: no changes in normalization
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_ref, consensus_alt = consensus.normalize("T", "C")
    assert consensus_ref == "T"
    assert consensus_alt == "C"


def test_normalization_case_3():
    '''
    test consensus.normalize() with identical ref/alt of length  < 1
    expectation: no ref, alt = first character
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_ref, consensus_alt = consensus.normalize("TAC", "TAC")
    assert consensus_ref == "T"
    assert consensus_alt == "T"


def test_normalization_case_4():
    '''
    test consensus.normalize() with differing ref/alt  with shared first character
    expectation: no changes in  normalization
    Parameters
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize("T", "TAC")
    assert consensus_ref == "T"
    assert consensus_alt == "TAC"


def test_normalization_case_5():
    '''
    test consensxus.normalize() with differing ref/alt  with shared first character
    expectation: normalized = shared first characters + all characters after 1st index
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_ref, consensus_alt = consensus.normalize("AT", "ATT")
    assert consensus_ref == "A"
    assert consensus_alt == "AT"


def test_normalization_case_6():
    '''
    test consensxus.normalize() GTA/GTATA
    expectation: G, GTA; first character matches, second 2 of ref matches to 
    alt, so alt becomes first  character + last 2
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize("GTA", "GTATA")
    assert consensus_ref == "G"
    assert consensus_alt == "GTA"


def test_normalization_case_7():
    '''
    test consensxus.normalize() GTA/GTATA
    expectation: G, GTA; first character matches, second 2 of ref matches to 
    alt, so alt becomes first  character + last 2
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize("GTA", "GTATA")
    assert consensus_ref == "G"
    assert consensus_alt == "GTA"



####################################
# test consensus.write_vcf() #
####################################

def test_write_vcf_case_1(testdir):
    '''
    test  consensus.write_vcf()
    Parameters
    testdir: str directory
    ---------
    Returns 
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    writeable_record = [chrom,pos,ref,alt] + record

    vcf = os.path.join(testdir, "vcf")
    counts = os.path.join(testdir, "counts")

    consensus.write_vcf([writeable_record], vcf, counts)

    assert os.path.exists(vcf)
    assert os.path.exists(counts)

    vcf = pd.read_csv(vcf, sep="\t")
    vcf = vcf.astype({"#CHROM": str, "FILTER":str})

    assert vcf.columns.tolist() == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    assert len(vcf) == 1

    assert vcf.values.tolist()[0] == [writeable_record[0], writeable_record[1], 
        writeable_record[4], writeable_record[2], writeable_record[3],
        writeable_record[5], writeable_record[6], '.']

    counts = pd.read_csv(counts, sep="\t")
    counts = counts.astype({"chrom": str})

    assert counts.columns.tolist() == ["chrom", "pos", "ID", "TR", "TA", "TD", "NR", "NA", "ND"]

    assert counts.values.tolist()[0]  == [writeable_record[0], writeable_record[1], writeable_record[4], 
    writeable_record[7], writeable_record[8], writeable_record[9], writeable_record[10],
                    writeable_record[11], writeable_record[12]]

    pypeliner.commandline.execute("rm", os.path.join(testdir, "vcf"))
    pypeliner.commandline.execute("rm", os.path.join(testdir, "counts"))




####################################
# test consensus.get_counts() #
####################################

def test_get_counts_case_1():
    '''
    test consensus.get_counts() with museq
    Parameters
    ---------
    Returns
    -------
    '''
    test1,_  = _get_test_model_call()
    assert consensus.get_counts(test1, "museq_snv", "normal", "tumor", test1.REF, test1.ALT) == (1, [2], 1, 1, [2], 1)


def test_get_counts_case_2():
    '''
    test consensus.get_counts() with museq
    Parameters
    ---------
    Returns
    -------
    '''
    test1,_  = _get_test_model_call()

    assert consensus.get_counts(test1, "strelka_snv", "normal", "tumor", test1.REF, test1.ALT) == (5, [6], 1, 5, [6], 1)


def test_get_counts_case_3():
    '''
    test consensus.get_counts() with rtg
    Parameters
    ---------
    Returns
    -------
    '''
    test1,_  = _get_test_model_call()

    assert consensus.get_counts(test1, "strelka_indel", "normal", "tumor", test1.REF, test1.ALT) == (3, [4], 1, 3, [4], 1)


def test_get_counts_case_4():
    '''
    test consensus.get_counts() with samtools
    Parameters
    ---------
    Returns
    -------
    '''
    test1,_  = _get_test_model_call()

    assert consensus.get_counts(test1, "mutect", "normal", "tumor", test1.REF, test1.ALT) == (1, [2], 1, 1, [2], 1)


