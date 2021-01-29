import random
import vcf
import gzip
import os
import consensus
import pypeliner
from shutil import copyfile
import csv
import pandas as pd
from collections import namedtuple


def _check_record(record, df):
    for i, row in df.iterrows():
        assert row[["chrom", "pos", "ref", "alt", "qual", "filter", "nr", "na", "nd", "id_counter"]].values.tolist() == record


def _get_test_record():
    '''
    make a dummy, correctly formatted snv with correct data
    Parameters
    ---------
    Returns
    -------
    '''
    alt = random.choice(["A", "T", "C", "G"])
    ref = random.choice(["A", "T", "C", "G"])
    chrom = random.choice(list(map(str, list(range(1, 23)) + ["X"])))
    pos = random.randint(0, 100000)

    test_record = [1, '.', 'nan', 1, 2, 3]
    return chrom, pos, ref, alt, test_record


def _get_test_model_call(DP = 1, RC=2,AC=2):
    #make example  call
    Call = namedtuple("Call", "DP RO AO RC AC AD")
    sample_data = vcf.model._Call(None, "sample_label", Call(DP,1,2,RC,AC, [1,2]) )

    freebayes_museq_rtg_example = vcf.model._Record("1", 10, 1, "A", [None], 1, 1, 
        None, None, None, [sample_data])

    samtools_example = vcf.model._Record("1", 10, 1, "A", [None], 1, 1, {"DP":1}, 
        None, None, [sample_data])

    return freebayes_museq_rtg_example, samtools_example

##################################
# test consensus.snv_consensus() #
##################################

def test_snv_consensus_case_1():
    '''
    test  consensus.snv_consensus() with identical record between museq/freebayes
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    freebayes = test_record
    rtg = []
    samtools = []

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})

    assert len(consensus_data) == 1

    _check_record([chrom, pos, ref, alt] + record, consensus_data )



def test_snv_consensus_case_2():
    '''
    test  consensus.snv_consensus() with identical record between museq/samtools
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    freebayes = []
    rtg = []
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_3():
    '''
    test  consensus.snv_consensus() with identical record between museq/rtg
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    freebayes = []
    rtg = test_record
    samtools = []

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_4():
    '''
    test  consensus.snv_consensus() with identical record between freebayes/samtools
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = []
    freebayes = test_record
    rtg = []
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_5():
    '''
    test  consensus.snv_consensus() with identical record between freebayes/rtg
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = []
    freebayes = test_record
    rtg = test_record
    samtools = []

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_6():
    '''
    test  consensus.snv_consensus() with identical record between samtools/freebayes
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = []
    freebayes = test_record
    rtg = []
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_7():
    '''
    test  consensus.snv_consensus() with identical record between samtools/rtg
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = []
    freebayes = []
    rtg = test_record
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


def test_snv_consensus_case_8():
    '''
    test  consensus.snv_consensus() with no identical records
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    museq = {(chrom, pos + 1, ref, alt): record}

    samtools = {(chrom, pos + 2, ref, alt): record}

    freebayes = {(chrom, pos + 3, ref, alt): record}

    rtg = {(chrom, pos + 4, ref, alt): record}


    consensus_data = consensus.snv_consensus(museq, samtools, freebayes, rtg)

    assert consensus_data == []


def test_snv_consensus_case_9():
    '''
    test  consensus.snv_consensus() with empty data
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_data = consensus.snv_consensus([], [], [], [])

    assert consensus_data == []


def test_snv_consensus_case_10():
    '''
    test  consensus.snv_consensus() with identical record between samtools/freebayes
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    freebayes = test_record
    rtg = []
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )



def test_snv_consensus_case_11():
    '''
    test  consensus.snv_consensus() with identical record between samtools/freebayes
    Parameters
    ---------
    Returns
    -------
    '''

    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    museq = test_record
    freebayes = test_record
    rtg = test_record
    samtools = test_record

    consensus_data = consensus.snv_consensus(museq, freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom": "str"})
    
    _check_record([chrom, pos, ref, alt] + record, consensus_data )


####################################
# test consensus.indel_consensus() #
####################################


def test_indel_consensus_case_1():
    '''
    test consensus.indel_consensus() with identical record between rtg/freebayes
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    freebayes = test_record
    rtg = test_record
    samtools = []

    consensus_data = consensus.indel_consensus(freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom":"str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_indel_consensus_case_2():
    '''
    test consensus.indel_consensus() with identical record between rtg/samtools
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    freebayes = []
    rtg = test_record
    samtools = test_record

    consensus_data = consensus.indel_consensus(freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom":"str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_indel_consensus_case_3():
    '''
    test consensus.indel_consensus() with identical record between freebayes/samtools
    Parameters
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    test_record = {(chrom, pos, ref, alt): record}

    freebayes = test_record
    rtg = []
    samtools = test_record

    consensus_data = consensus.indel_consensus(freebayes, rtg, samtools)
    
    consensus_data = pd.DataFrame(consensus_data, 
        columns = ["chrom", "pos", "ref", "alt", "id_counter", "qual", "filter", "nr", "na", "nd"]
    )

    consensus_data = consensus_data.astype({"chrom":"str"})
    
    assert not consensus_data[(consensus_data.pos==pos) & (consensus_data["chrom"]==chrom)].empty


def test_indel_consensus_case_4():
    '''
    test  consensus.snv_consensus() with empty data
    Parameters
    testdir: str directory
    ---------
    Returns
    -------
    '''
    consensus_data = consensus.indel_consensus([], [], [])

    assert consensus_data == []


def test_indel_consensus_case_5():
    '''
    test  consensus.snv_consensus() with no shared records
    Parameters
    testdir: str directory
    ---------
    Returns
    -------
    '''
    chrom, pos, ref, alt, record = _get_test_record()

    samtools = {(chrom, pos + 2, ref, alt): record}

    freebayes = {(chrom, pos + 3, ref, alt): record}

    rtg = {(chrom, pos + 4, ref, alt): record}

    consensus_data = consensus.indel_consensus(samtools, freebayes, rtg)

    assert consensus_data == []


####################################
# test consensus.normalize() #
####################################


def test_normalization_case_1():
    '''
    test consensus.normalize() with None data
    expectation: None outputs
    Parameters
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize(None, None)
    assert consensus_ref == None
    assert consensus_alt == None


def test_normalization_case_2():
    '''
    test consensus.normalize() with empty data
    expectation: empty outputs
    Parameters
    ---------
    Returns
    -------
    '''

    consensus_ref, consensus_alt = consensus.normalize("", "")
    assert consensus_ref == ""
    assert consensus_alt == ""


def test_normalization_case_3():
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


def test_normalization_case_4():
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


def test_normalization_case_5():
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


def test_normalization_case_6():
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


def test_normalization_case_7():
    '''
    test consensxus.normalize() with differing ref/alt  with shared first character
    expectation: normalized = shared first characters + all characters after 1st index
    Parameters
    ---------
    Returns
    -------
    '''
    consensus_ref, consensus_alt = consensus.normalize("TGACCAT", "TACTCAA")
    assert consensus_ref == "TACCAT"
    assert consensus_alt == "TCTCAA"



def test_normalization_case_8():
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


def test_normalization_case_9():
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

    assert counts.columns.tolist() == ["chrom", "pos", "ID", "NR", "NA", "ND"]

    assert counts.values.tolist()[0]  == [writeable_record[0], writeable_record[1], 
        writeable_record[4], writeable_record[7], writeable_record[8],
        writeable_record[9]]

    pypeliner.commandline.execute("rm", os.path.join(testdir, "vcf"))
    pypeliner.commandline.execute("rm", os.path.join(testdir, "counts"))


####################################
# test consensus.get_counts() #
####################################

def test_get_counts_case_1():
    '''
    test consensus.get_counts() with freebayes
    Parameters
    ---------
    Returns
    -------
    '''
    freebayes_museq_rtg_example, _ = _get_test_model_call()

    assert consensus.get_counts(freebayes_museq_rtg_example, "freebayes", "sample_label") == (1, [2], 1)

def test_get_counts_case_2():
    '''
    test consensus.get_counts() with museq
    Parameters
    ---------
    Returns
    -------
    '''
    freebayes_museq_rtg_example, _ = _get_test_model_call()
    
    assert consensus.get_counts(freebayes_museq_rtg_example, "museq_germline", "sample_label") == (1, [2], 1)

def test_get_counts_case_3():
    '''
    test consensus.get_counts() with rtg
    Parameters
    ---------
    Returns
    -------
    '''
    freebayes_museq_rtg_example, _ = _get_test_model_call()
    
    assert consensus.get_counts(freebayes_museq_rtg_example, "rtg", "sample_label") == (1, [2], 1)

def test_get_counts_case_4():
    '''
    test consensus.get_counts() with samtools
    Parameters
    ---------
    Returns
    -------
    '''
    _, samtools_example = _get_test_model_call()
    
    assert consensus.get_counts(samtools_example, "samtools", "sample_label") == ("NA", ["NA"], 1)

test_normalization_case_7()


