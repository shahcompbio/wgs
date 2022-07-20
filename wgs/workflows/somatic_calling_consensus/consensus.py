from collections import defaultdict

import vcf


def get_reader(filename):
    return vcf.Reader(filename=filename)


def get_counts(record, caller, tumor_id, normal_id, ref, alts):
    normal = [v for v in record.samples if v.sample == normal_id]
    assert len(normal) == 1
    normal = normal[0]

    tumor = [v for v in record.samples if v.sample == tumor_id]
    assert len(tumor) == 1
    tumor = tumor[0]

    if caller == 'museq_snv':
        tumor_depth = tumor['DP']
        tumor_ref = tumor['RC']
        tumor_alt = [tumor['AC']]
        normal_depth = normal['DP']
        normal_ref = normal['RC']
        normal_alt = [normal['AC']]
    elif caller == 'strelka_snv':
        assert len(alts) == 1
        alt = str(alts[0])
        tumor_depth = tumor['DP']
        tumor_ref = tumor[ref + 'U'][0]
        tumor_alt = [tumor[alt + 'U'][0]]
        normal_depth = normal['DP']
        normal_ref = normal[ref + 'U'][0]
        normal_alt = [normal[alt + 'U'][0]]
    elif caller == 'strelka_indel':
        assert len(alts) == 1
        tumor_depth = tumor['DP']
        tumor_ref = tumor['TAR'][0]
        tumor_alt = [tumor['TIR'][0]]
        normal_depth = normal['DP']
        normal_ref = normal['TAR'][0]
        normal_alt = [normal['TIR'][0]]
    elif caller == 'mutect':
        tumor_depth = tumor['DP']
        tumor_ref = tumor['AD'][0]
        tumor_alt = tumor['AD'][1:]
        normal_depth = normal['DP']
        normal_ref = normal['AD'][0]
        normal_alt = normal['AD'][1:]

    else:
        raise NotImplementedError()

    return tumor_ref, tumor_alt, tumor_depth, normal_ref, normal_alt, normal_depth


def fetch_vcf(filename, chromosome, caller):
    snv_data = {}
    indel_data = {}

    vcf_reader = get_reader(filename)

    tumor_sample = vcf_reader.metadata['tumor_sample'][0]
    normal_sample = vcf_reader.metadata['normal_sample'][0]

    id_counter = 0

    try:
        records = vcf_reader.fetch(chromosome)
    except ValueError:
        return snv_data, indel_data

    for record in records:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alts = record.ALT
        vcf_filter = record.FILTER

        if caller == 'mutect' and vcf_filter:
            continue
        elif vcf_filter is None:
            vcf_filter = '.'
        elif vcf_filter == []:
            vcf_filter = 'PASS'
        else:
            assert len(vcf_filter) <= 1
            vcf_filter = vcf_filter[0]

        if alts == [None]:
            continue
        tr, tas, td, nr, nas, nd = get_counts(record, caller, tumor_sample, normal_sample, ref, alts)

        assert len(alts) == len(tas) == len(nas)

        for (alt, ta, na) in zip(alts, tas, nas):
            alt = str(alt)

            data = [record.QUAL, vcf_filter, tr, ta, td, nr, na, nd, '{}_{}'.format(caller, id_counter)]

            if len(ref) == len(alt):
                for i, (rb, ab) in enumerate(zip(ref, alt)):
                    if not rb == ab:
                        snv_data[(chrom, pos + i, rb, ab)] = data
                        id_counter += 1
            else:
                indel_data[(chrom, pos)] = (data, ref, alt)
                id_counter += 1

    return snv_data, indel_data


def snv_consensus(museq, strelka, mutect):
    outdata = defaultdict(int)
    for v in museq:
        outdata[v] += 1
    for v in strelka:
        outdata[v] += 1
    for v in mutect:
        outdata[v] += 1

    consensus = []
    for k in outdata:
        if outdata[k] <= 1:
            continue

        if k in museq:
            qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = museq[k]
        else:
            qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = strelka[k]

        consensus.append([
            k[0], k[1], k[2], k[3], id_count, qual, vcf_filter, tr, ta, td, nr, na, nd
        ])

    return consensus


def indel_consensus(strelka_indel, mutect_indel):
    consensus = []

    for k in strelka_indel:
        if k in mutect_indel:
            chrom, pos = k
            mutect_data, mutect_ref, mutect_alt = mutect_indel[k]
            strelka_data, strelka_ref, strelka_alt = strelka_indel[k]

            if mutect_ref == strelka_ref and mutect_alt == strelka_alt:
                qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = mutect_data
                consensus.append(
                    [chrom, pos, mutect_ref, mutect_alt, id_count, qual, vcf_filter, tr, ta, td, nr, na, nd])
            else:
                qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = mutect_data
                consensus.append(
                    [chrom, pos, mutect_ref, mutect_alt, id_count, qual, vcf_filter, tr, ta, td, nr, na, nd])
                qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = strelka_data
                consensus.append(
                    [chrom, pos, strelka_ref, strelka_alt, id_count, qual, vcf_filter, tr, ta, td, nr, na, nd])
        else:
            chrom, pos = k
            strelka_data, strelka_ref, strelka_alt = strelka_indel[k]
            qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = strelka_data
            consensus.append([chrom, pos, strelka_ref, strelka_alt, id_count, qual, vcf_filter, tr, ta, td, nr, na, nd])

    for k in mutect_indel:
        if k not in strelka_indel:
            chrom, pos = k
            mutect_data, mutect_ref, mutect_alt = mutect_indel[k]
            qual, vcf_filter, tr, ta, td, nr, na, nd, id_count = mutect_data
            consensus.append([chrom, pos, mutect_ref, mutect_alt, id_count, qual, vcf_filter, tr, ta, td, nr, na, nd])

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
        museq_snv_vcf,
        strelka_snv_vcf,
        mutect_snv_vcf,
        strelka_indel_vcf,
        consensus_vcf,
        counts_output,
        chromosomes,
):
    with open(consensus_vcf, 'wt') as writer:
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    with open(counts_output, 'wt') as writer:
        writer.write("chrom\tpos\tID\tTR\tTA\tTD\tNR\tNA\tND\n")

    for chromosome in chromosomes:
        museq_calls, _ = fetch_vcf(museq_snv_vcf, chromosome, 'museq_snv')
        strelka_snv, _ = fetch_vcf(strelka_snv_vcf, chromosome, 'strelka_snv')
        _, strelka_indel = fetch_vcf(strelka_indel_vcf, chromosome, 'strelka_indel')
        mutect_snv, mutect_indel = fetch_vcf(mutect_snv_vcf, chromosome, 'mutect')

        consensus = snv_consensus(museq_calls, strelka_snv, mutect_snv)
        consensus += indel_consensus(strelka_indel, mutect_indel)

        write_vcf(consensus, consensus_vcf, counts_output)
