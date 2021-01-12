from collections import defaultdict

import vcf


def get_reader(filename):
    return vcf.Reader(filename=filename)


def get_counts(record, caller, sample_id):
    sample = [v for v in record.samples if v.sample == sample_id]
    assert len(sample) == 1, (sample, caller, sample_id, record, sample_id)
    sample = sample[0]

    if caller == 'museq_germline':
        depth = sample['DP']
        ref = sample['RC']
        alt = sample['AC']
    elif caller == 'freebayes':
        depth = sample['DP']
        ref = sample['RO']
        alt = sample['AO']
    elif caller == 'rtg':
        depth = sample['DP']
        ad = sample['AD']
        if isinstance(ad, list):
            assert len(ad) > 1
            ref = ad[0]
            alt = ad[1:]
        else:
            assert record.ALT == [None]
            assert isinstance(ad, int)
            ref = ad
            alt = ['NA']
            raise Exception('TODO')
    elif caller == 'samtools':
        depth = record.INFO['DP']
        ref = 'NA'
        alt = ['NA']
    else:
        raise NotImplementedError()
    if isinstance(alt, int):
        alt = [alt]

    return ref, alt, depth


def fetch_vcf(filename, chromosome, caller):
    snv_data = {}
    indel_data = {}

    vcf_reader = get_reader(filename)

    sample_id = vcf_reader.metadata['normal_sample'][0]

    try:
        records = vcf_reader.fetch(chromosome)
    except ValueError:
        return snv_data, indel_data

    id_counter = 0

    for record in records:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alts = record.ALT
        filter = record.FILTER

        if not filter:
            filter = '.'
        else:
            assert len(filter) <= 1
            filter = filter[0]

        ref_count, alt_counts, depth = get_counts(record, caller, sample_id)

        for alt, alt_count in zip(alts, alt_counts):
            alt = str(alt)
            data = [record.QUAL, filter, ref_count, alt_count, depth, '{}_{}'.format(caller, id_counter)]
            if len(ref) == len(alt):
                for i, (rb, ab) in enumerate(zip(ref, alt)):
                    if not rb == ab:
                        snv_data[(chrom, pos + i, rb, ab)] = data
                        id_counter += 1
            else:
                ref, alt = normalize(ref, alt)
                indel_data[(chrom, pos, ref, alt)] = data
                id_counter += 1

    return snv_data, indel_data


def snv_consensus(museq, freebayes, rtg, samtools):
    outdata = defaultdict(int)
    for v in museq:
        outdata[v] += 1
    for v in freebayes:
        outdata[v] += 1
    for v in rtg:
        outdata[v] += 1
    for v in samtools:
        outdata[v] += 1

    consensus = []
    for k in outdata:
        if outdata[k] <= 1:
            continue

        if k in museq:
            qual, filter, nr, na, nd, id_counter = museq[k]
        elif k in freebayes:
            qual, filter, nr, na, nd, id_counter = freebayes[k]
        else:
            qual, filter, nr, na, nd, id_counter = rtg[k]

        consensus.append([
            k[0], k[1], k[2], k[3], id_counter, qual, filter, nr, na, nd
        ])

    return consensus


def normalize(ref, alt):
    if not ref and not alt:
        return ref, alt

    if len(ref) == 1 or len(alt) == 1:
        return ref, alt

    if not ref[0] == alt[0]:
        return ref, alt

    common_base = ref[0]
    ref = ref[1:]
    alt = alt[1:]

    for i, (ref_v, alt_v) in enumerate(zip(ref[::-1], alt[::-1])):
        if not ref_v == alt_v:
            break

    ref = ref[i + 1:]
    alt = alt[i + 1:]

    ref = common_base + ref
    alt = common_base + alt

    return ref, alt


def indel_consensus(freebayes, rtg, samtools):
    outdata = defaultdict(int)
    for v in freebayes:
        outdata[v] += 1
    for v in rtg:
        outdata[v] += 1
    for v in samtools:
        outdata[v] += 1

    consensus = []

    for k in outdata:
        if outdata[k] <= 1:
            continue
        chrom, pos, ref, alt = k

        fb_data = freebayes.get(k)
        rtg_data = rtg.get(k)

        if fb_data:
            qual, vcf_filter, nr, na, nd, id_counter = fb_data
        else:
            qual, vcf_filter, nr, na, nd, id_counter = rtg_data

        consensus.append((chrom, pos, ref, alt, id_counter, qual, vcf_filter, nr, na, nd))

    return consensus


def write_vcf(consensus, vcf_output, counts_output):
    with open(vcf_output, 'a') as outfile, open(counts_output, 'a') as count_file:
        for call in consensus:
            outstr = [call[0], str(call[1]), str(call[4]), call[2], call[3], str(call[5]), call[6], '.']
            outstr = '\t'.join(outstr) + '\n'
            outfile.write(outstr)
            outstr = [call[0], str(call[1]), str(call[4]), str(call[7]), str(call[8]), str(call[9])]
            outstr = '\t'.join(outstr) + '\n'
            count_file.write(outstr)


def main(
        museq_vcf,
        freebayes_vcf,
        rtg_vcf,
        samtools_vcf,
        consensus_vcf,
        counts_output,
        chromosomes
):
    for chromosome in chromosomes:
        museq_calls, _ = fetch_vcf(museq_vcf, chromosome, 'museq_germline')
        freebayes_snv, freebayes_indel = fetch_vcf(freebayes_vcf, chromosome, 'freebayes')
        rtg_snv, rtg_indel = fetch_vcf(rtg_vcf, chromosome, 'rtg')
        samtools_snv, samtools_indel = fetch_vcf(samtools_vcf, chromosome, 'samtools')

        consensus = snv_consensus(museq_calls, freebayes_snv, rtg_snv, samtools_snv)
        consensus += indel_consensus(freebayes_indel, rtg_indel, samtools_indel)

        write_vcf(consensus, consensus_vcf, counts_output)
