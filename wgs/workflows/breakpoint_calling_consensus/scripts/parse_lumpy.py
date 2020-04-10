import pandas as pd
import vcf

VCF_FILE = "lumpy.vcf"


def get_reader(vcf_file):
    return vcf.Reader(open(vcf_file, 'rt'))


def parse_vcf(vcf_file):
    reader = get_reader(vcf_file)

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
            if isinstance(v, list):
                v = ';'.join(map(str, v))
            data[k] = v

        for sample in record.samples:
            sample_name = sample.sample
            sample_data = sample.data
            for k, v in sample_data._asdict().items():
                if isinstance(v, list):
                    v = ';'.join([str(val) for val in v])
                k = '{}_{}'.format(sample_name, k)
                data[k] = v

        yield data


def parse_vcf_group(data):
    out_record = []

    for record in data:
        if record['SVTYPE'] == 'BND':
            out_record.append(record)
            if len(out_record) == 2:
                yield out_record
                out_record = []
        else:
            yield (record,)


def create_data(vcfdata):
    same_cols = [
        'ref', 'qual', 'SVTYPE', 'SVLEN',
        'IMPRECISE', 'SU', 'PE', 'SR', 'GT',
        'SU', 'PE', 'SR'
    ]

    dup_cols = [
        'alt', 'STRANDS', 'CIPOS', 'CIEND', 'CIPOS95', 'CIEND95',
    ]

    for record in vcfdata:

        out_record = {}

        if len(record) == 1:
            record = record[0]

            out_record['chromosome_1'] = record['chrom']
            out_record['chromosome_2'] = record['chrom']

            out_record['position_1'] = record['pos']
            out_record['position_2'] = record['END']

            strands = record['STRANDS'].split(':')[0]
            assert len(strands) == 2
            out_record['strand_1'] = strands[0]
            out_record['strand_2'] = strands[1]

            for col in dup_cols:
                value = record[col] if col in record else None
                out_record[col + '_1'] = value
                out_record[col + '_2'] = value

            for col in same_cols:
                value = record[col] if col in record else None
                out_record[col] = value

        elif len(record) == 2:
            if record[0]['MATEID'].endswith('_1'):
                mate1 = record[0]
                mate2 = record[1]
            else:
                mate1 = record[1]
                mate2 = record[0]

            out_record['chromosome_1'] = mate1['chrom']
            out_record['chromosome_2'] = mate2['chrom']

            out_record['position_1'] = mate1['pos']
            out_record['position_2'] = mate2['pos']

            strands_1 = mate1['STRANDS'].split(':')[0]
            strands_2 = mate2['STRANDS'].split(':')[0]
            assert strands_1 == strands_2[::-1]
            out_record['strand_1'] = strands_1[0]
            out_record['strand_2'] = strands_2[0]

            for col in dup_cols:
                value_1 = mate1[col] if col in record[0] else None
                value_2 = mate2[col] if col in record[1] else None
                out_record[col + '_1'] = value_1
                out_record[col + '_2'] = value_2

            for col in same_cols:
                value = mate1[col] if col in mate1 else None
                out_record[col] = value

                if value:
                    assert mate1[col] == mate2[col]
        else:
            raise Exception('unknown record')

        yield out_record


def convert_to_df(data):
    data = list(data)
    df = pd.DataFrame(data)

    return df


def filter_calls(data, filters):
    if filters.get('tumour_read_support_threshold'):
        data = data[data['SU'] >= filters['tumour_read_support_threshold']]

    if filters.get('deletion_size_threshold'):
        data = data[(data['SVTYPE'] == 'DEL') & (data['SVLEN'] >= filters['deletion_size_threshold'])]

    if filters.get('chromosomes'):
        data = data[data['chromosome_1'].isin(filters['chromosomes'])]
        data = data[data['chromosome_2'].isin(filters['chromosomes'])]

    return data


def write(output, data):
    data.to_csv(output, index=False)


