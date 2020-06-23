import vcf


def read_ref_positions(positions):
    ref_pos = set()
    freader = open(positions)
    for line in freader:
        line = line.strip().split(':')

        chrom = line[0]
        pos = int(line[1])

        ref_pos.add((chrom, pos))

    return ref_pos


def get_reader(filename):
    return vcf.Reader(filename=filename)


def vcf_to_counts(filename, outfile, ref_positions):
    vcf_reader = get_reader(filename)

    ref_positions = read_ref_positions(ref_positions)

    tumor_sample = vcf_reader.metadata['tumor_sample'][0]

    with open(outfile, 'wt') as writer:

        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF

            if (chrom, pos) not in ref_positions:
                continue

            tumor = [v for v in record.samples if v.sample == tumor_sample]
            assert len(tumor) == 1
            tumor = tumor[0]

            tumor_ref = str(tumor['RC'])
            tumor_alt = str(tumor['AC'])
            pos = str(pos)

            outstr = '\t'.join([chrom, pos, ref, tumor_ref, 'X', tumor_alt]) + '\n'

            writer.write(outstr)
