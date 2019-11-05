import vcf
from wgs.utils import helpers

VCF_FILE = "museq_paired_annotated.vcf"

ANNOTATIONS = ['ann', 'ma', 'dbsnp', 'cosmic', 'lof', 'nmd']


def get_reader(vcf_file):
    return vcf.Reader(open(vcf_file, 'rt'))


def get_snpeff_cols(reader):
    info = reader.infos['ANN']

    desc = info.desc

    desc = desc.split("'")[1]

    desc = desc.split('|')

    desc = [v.replace(' ', '') for v in desc]

    return desc


def parse_snpeff(snpeff_cols, snpeff_entries, chrom, pos):
    records = []
    for record in snpeff_entries:
        record = record.strip().split('|')

        record = dict(zip(snpeff_cols, record))

        record['chrom'] = chrom
        record['pos'] = pos

        records.append(record)

    return records


def parse_vcf(vcf_file):
    reader = get_reader(vcf_file)

    snpeff_cols = get_snpeff_cols(reader)

    for record in reader:
        data = {
            'chrom': record.CHROM,
            'pos': record.POS,
            'ref': record.REF,
            'alt': record.ALT,
            'qual': record.QUAL,
        }

        info = record.INFO

        for k, v in info.items():
            if k.lower() in ANNOTATIONS:
                continue
            if isinstance(v, list):
                v = ';'.join(v)
            data[k] = v

        for sample in record.samples:

            sample_type = sample.sample
            sample_data = sample.data

            for k, v in sample_data._asdict().items():
                k += "_" + sample_type
                if isinstance(v, list):
                    v = ';'.join([str(val) for val in v])
                data[k] = v

        snpeff_annotations = parse_snpeff(snpeff_cols, info['ANN'], record.CHROM, record.POS)
        # annotations['dbsnp'] = info['DBSNP']
        # annotations['ma'] = info['MA']
        # annotations['cosmic'] = info['Cosmic']

        yield data, (snpeff_annotations,)


def write_record(record, header, outfile):
    outstr = [record[col] for col in header]

    outstr = ','.join(map(str, outstr)) + '\n'

    outfile.write(outstr)


def write(outputs, data):
    with helpers.GetFileHandle(outputs['primary'], 'w') as primary_outfile, \
            helpers.GetFileHandle(outputs['snpeff'], 'w') as snpeff_outfile:

        primary_header = None
        snpeff_header = None

        for record in data:
            primary, annotations = record

            if not primary_header:
                primary_header = list(primary.keys())

            write_record(primary, primary_header, primary_outfile)

            snpeff_anns = annotations[0]

            for snpeff_record in snpeff_anns:
                if not snpeff_header:
                    snpeff_header = list(snpeff_record.keys())

                write_record(snpeff_record, snpeff_header, snpeff_outfile)


if __name__ == "__main__":
    vcfdata = parse_vcf(VCF_FILE)

    write({'primary': 'parsed', 'snpeff': 'snpeff'}, vcfdata)
