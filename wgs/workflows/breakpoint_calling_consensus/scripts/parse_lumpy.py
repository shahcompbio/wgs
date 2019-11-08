import pandas as pd
import vcf

VCF_FILE = "lumpy.vcf"

ANNOTATIONS = []


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
            if k.lower() in ANNOTATIONS:
                continue
            if isinstance(v, list):
                v = ';'.join(map(str, v))
            data[k] = v

        for sample in record.samples:

            sample_type = sample.sample
            sample_data = sample.data

            for k, v in sample_data._asdict().items():
                k += "_" + sample_type
                if isinstance(v, list):
                    v = ';'.join([str(val) for val in v])
                data[k] = v

        yield data


def convert_to_df(data):
    data = list(data)
    df = pd.DataFrame(data)

    return df


def write(output, data):
    df = convert_to_df(data)
    df.to_csv(output, index=False)


if __name__ == "__main__":
    vcfdata = parse_vcf(VCF_FILE)

    write('parsed.csv', vcfdata)
