import os
from collections import defaultdict

import vcf


def get_reader(filename):
    return vcf.Reader(filename=filename)


def get_counts(record, caller, tumor_id, normal_id, ref, alts):
    if caller == 'museq':
        normal = [v for v in record.samples if v.sample == normal_id]
        assert len(normal) == 1
        normal = normal[0]
        normal_depth = normal['DP']
        normal_ref = normal['RC']
        normal_alt = [normal['AC']]
    elif caller == 'samtools':
        normal_depth = record.INFO['DP']
        normal_ref = 'NA'
        normal_alt = ['NA' for alt in alts]
    else:
        raise NotImplementedError()

    return normal_ref, normal_alt, normal_depth


def fetch_vcf(filename, chromosome, caller):
    snv_data = {}
    indel_data = {}

    vcf_reader = get_reader(filename)

    id_counter = 0

    try:
        records = vcf_reader.fetch(chromosome)
    except ValueError:
        return snv_data, indel_data

    for record in records:
        chrom = record.CHROM

        assert chrom == chromosome
        pos = record.POS
        ref = record.REF
        alts = record.ALT
        filter_call = record.FILTER

        if caller == 'museq' and record.INFO['PR'] < 0.8:
            continue

        if not filter_call:
            filter_call = '.'
        else:
            assert len(filter_call) <= 1
            filter_call = filter_call[0]

        if caller == 'museq' and record.INFO < 0.85:
            continue

        tum_id = 'TUMOR' if caller == 'samtools' else 'TUMOUR'

        nr, nas, nd = get_counts(record, caller, tum_id, 'NORMAL', ref, alts)

        assert len(alts) == len(nas)

        for (alt, na) in zip(alts, nas):
            alt = str(alt)

            data = [record.QUAL, filter_call, 'NA', 'NA', 'NA', nr, na, nd, "{}_{}".format(caller, id_counter)]

            if len(ref) == len(alt):
                for i, (rb, ab) in enumerate(zip(ref, alt)):
                    if not rb == ab:
                        snv_data[(chrom, pos + i, rb, ab)] = data
                        id_counter += 1
            else:
                indel_data[(chrom, pos)] = (data, ref, alt)
                id_counter += 1

    return snv_data, indel_data


def snv_consensus(museq, samtools):
    outdata = defaultdict(int)
    for v in museq:
        outdata[v] += 1
    for v in samtools:
        outdata[v] += 1

    consensus = []
    for k in outdata:
        if outdata[k] <= 1:
            continue

        if k in museq:
            qual, filter_call, tr, ta, td, nr, na, nd, id_count = museq[k]
        else:
            qual, filter_call, tr, ta, td, nr, na, nd, id_count = samtools[k]

        consensus.append([
            k[0], k[1], k[2], k[3], id_count, qual, filter_call, tr, ta, td, nr, na, nd
        ])

    return consensus


def normalize(ref, alt):
    if len(ref) == 1 or len(alt) == 1:
        return ref, alt

    assert ref[0] == alt[0]
    common_base = ref[0]
    ref = ref[1:]
    alt = alt[1:]

    for i, (ref_v, alt_v) in enumerate(zip(ref[::-1], alt[::-1])):
        assert ref_v == alt_v

    ref = ref[i + 1:]
    alt = alt[i + 1:]

    ref = common_base + ref
    alt = common_base + alt

    return ref, alt


def indel_consensus(samtools_indel):
    consensus = []
    for k in samtools_indel:
        chrom, pos = k
        samtools_data, samtools_ref, samtools_alt = samtools_indel[k]
        samtools_ref, samtools_alt = normalize(samtools_ref, samtools_alt)
        qual, flter, tr, ta, td, nr, na, nd, id_count = samtools_data
        consensus.append([chrom, pos, samtools_ref, samtools_alt, id_count, qual, flter, tr, ta, td, nr, na, nd])
    return consensus


def write_vcf(consensus, vcf_output, counts_output):
    with open(vcf_output, 'a') as outfile, open(counts_output, 'a') as count_file:
        for call in consensus:
            outstr = [call[0], str(call[1]), str(call[4]), call[2], call[3], str(call[5]), call[6], '.']
            outstr = '\t'.join(outstr) + '\n'
            outfile.write(outstr)

            outstr = [call[0], str(call[1]), str(call[4]), str(call[7]), str(call[8]), str(call[9]), str(call[10]),
                      str(call[11]), str(call[12])]
            outstr = '\t'.join(outstr) + '\n'
            if '[' in outstr:
                raise Exception(call)
            count_file.write(outstr)


def main(
        museq_vcf,
        samtools_vcf,
        consensus_vcf,
        counts_output,
        chromosomes,
):
    if os.path.exists(consensus_vcf):
        os.remove(consensus_vcf)
    if os.path.exists(counts_output):
        os.remove(counts_output)

    for chromosome in chromosomes:
        museq_calls, _ = fetch_vcf(museq_vcf, chromosome, 'museq')
        samtools_calls, samtools_indels = fetch_vcf(samtools_vcf, chromosome, 'samtools')

        consensus = snv_consensus(museq_calls, samtools_calls)
        consensus += indel_consensus(samtools_indels)

        write_vcf(consensus, consensus_vcf, counts_output)
