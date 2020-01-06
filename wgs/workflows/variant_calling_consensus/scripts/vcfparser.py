import vcf
from wgs.utils import helpers

VCF_FILE = "museq_single_annotated.vcf.gz"

ANNOTATIONS = ['ann', 'ma', 'dbsnp', 'cosmic', 'lof', 'nmd', '1000gen', 'low_mappability']


class VcfParser(object):
    def __init__(self, vcf_file, outfile, snpeff_outfile, ma_outfile, ids_outfile, filters):
        '''
        constructor for parser
        note, if filter_low_mappability is true,
        will look for a "fxblacklist" in the parser config
        '''
        self.vcf_file = vcf_file
        self.outfile = outfile
        self.snpeff_outfile = snpeff_outfile
        self.ma_outfile = ma_outfile
        self.ids_outfile = ids_outfile

        self.reader = self.get_reader(self.vcf_file)

        self.snpeff_cols, self.ma_cols, self.ids_cols, self.primary_cols = self.init_headers()

        self.cols = None

        self.filters = filters

    def __enter__(self):
        self.outfile = helpers.GetFileHandle(self.outfile, 'wt').handler
        self.snpeff_outfile = helpers.GetFileHandle(self.snpeff_outfile, 'wt').handler
        self.ma_outfile = helpers.GetFileHandle(self.ma_outfile, 'wt').handler
        self.ids_outfile = helpers.GetFileHandle(self.ids_outfile, 'wt').handler
        self.write_headers()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.outfile.close()
        self.snpeff_outfile.close()
        self.ma_outfile.close()
        self.ids_outfile.close()

    def init_headers(self):

        snpeff_cols = self.get_cols_from_header(self.reader, 'ANN') + ['chrom', 'pos']

        ma_cols = self.get_cols_from_header(self.reader, 'MA') + ['chrom', 'pos']

        primary_cols = self.get_primary_cols_from_header(self.reader)

        ids_cols = ['chrom', 'pos', 'value', 'type']

        return snpeff_cols, ma_cols, ids_cols, primary_cols

    def write_headers(self):
        self.write_header(self.snpeff_cols, self.snpeff_outfile)
        self.write_header(self.ma_cols, self.ma_outfile)
        self.write_header(self.ids_cols, self.ids_outfile)
        self.write_header(self.primary_cols, self.outfile)

    def get_reader(self, vcf_file):
        return vcf.Reader(filename=vcf_file)

    def get_primary_cols_from_header(self, reader):
        try:
            row1 = next(reader)
        except StopIteration:
            return ["chrom", "pos", "ref", "alt", "qual", "filter"]

        return list(self.parse_main_cols(row1).keys())

    def get_cols_from_header(self, reader, key):
        try:
            info = reader.infos[key]
        except KeyError:
            return ['chrom', 'pos']

        desc = info.desc

        desc = desc.split("'")
        if len(desc) == 1:
            desc = desc[0].split("(")[1]
        else:
            desc = desc[1]

        desc = desc.split('|')

        desc = [v.replace(' ', '') for v in desc]

        return desc

    def parse_snpeff(self, snpeff_cols, snpeff_entries, chrom, pos):
        records = []
        for record in snpeff_entries:
            record = record.strip().split('|')

            record = dict(zip(snpeff_cols, record))

            record['chrom'] = chrom
            record['pos'] = pos

            records.append(record)

        return records

    def parse_mutation_assessor(self, cols, record, chrom, pos):
        if not record:
            return
        record = record.strip().split('|')

        record = dict(zip(cols, record))

        record['chrom'] = chrom
        record['pos'] = pos

        return record

    def parse_list_annotations(self, records, chrom, pos, label):
        if records == [None]:
            return []

        parsed_records = []
        for entry in records:
            parsed_records.append(
                {'chrom': chrom, 'pos': pos, 'value': entry, 'type': label}
            )

        return parsed_records

    def parse_flag_annotation(self, record, chrom, pos, label):
        if not record:
            return []

        return [{'chrom': chrom, 'pos': pos, 'value': record, 'type': label}]

    def parse_main_cols(self, record):
        data = {
            'chrom': record.CHROM,
            'pos': record.POS,
            'ref': record.REF,
            'alt': record.ALT,
            'qual': record.QUAL,
        }

        data['alt'] = ';'.join(map(str, data['alt']))

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

        return data

    def eval_expr(self, val, operation, threshold):
        if operation == "gt":
            if val > threshold:
                return True
        elif operation == 'ge':
            if val >= threshold:
                return True
        elif operation == 'lt':
            if val < threshold:
                return True
        elif operation == 'le':
            if val <= threshold:
                return True
        elif operation == 'eq':
            if val == threshold:
                return True
        elif operation == 'ne':
            if not val == threshold:
                return True
        elif operation == 'in':
            if val in threshold:
                return True
        elif operation == 'notin':
            if not val in threshold:
                return True
        else:
            raise Exception("unknown operator type: {}".format(operation))

        return False

    def filter_records(self, record, annotations):

        for vcf_filter in self.filters:
            filter_name, relationship, value = vcf_filter

            if filter_name in record:
                if self.eval_expr(record[filter_name], relationship, value):
                    return True

            for annotation in annotations:
                if filter_name == annotation['type']:
                    if self.eval_expr(record[filter_name], relationship, value):
                        return True

    def parse_record(self, record):
        data = self.parse_main_cols(record)
        info = record.INFO

        snpeff_annotations = self.parse_snpeff(self.snpeff_cols, info['ANN'], record.CHROM, record.POS)
        ma_annotations = self.parse_mutation_assessor(self.ma_cols, info['MA'], record.CHROM, record.POS)

        id_annotations = self.parse_list_annotations(info['DBSNP'], record.CHROM, record.POS, 'dbsnp')
        id_annotations += self.parse_list_annotations(info['Cosmic'], record.CHROM, record.POS, 'cosmic')
        ## TODO: I expected to see a 1000gen: False in the record. but the key is missing.
        id_annotations += self.parse_flag_annotation(info.get('1000Gen'), record.CHROM, record.POS, '1000Gen')
        id_annotations += self.parse_flag_annotation(info.get('LOW_MAPPABILITY'), record.CHROM, record.POS,
                                                     'LOW_MAPPABILITY')
        return data, snpeff_annotations, id_annotations, ma_annotations

    def parse_vcf(self):
        for record in self.reader:

            data, snpeff_annotations, id_annotations, ma_annotations = self.parse_record(record)

            if self.filter_records(data, id_annotations):
                continue

            yield data, (snpeff_annotations, ma_annotations, id_annotations)

    def write_record(self, record, header, outfile):
        '''
        writes a record to outfile
        :param record: record to write
        :param header: header containg colnames
        :param outfile: path to outfile
        '''
        outstr = [record[col] for col in header]
        outstr = ','.join(map(str, outstr)) + '\n'
        outfile.write(outstr)

    def write_header(self, header, outfile):

        outstr = ','.join(map(str, header)) + '\n'
        outfile.write(outstr)

    def write(self):
        '''
        writes the parser to a csv
        :param blacklist:
        '''
        data = self.parse_vcf()
        primary_header = self.primary_cols
        for record in data:
            primary, annotations = record
            if not primary_header:
                primary_header = list(primary.keys())
                self.write_header(primary_header, self.outfile)
            self.write_record(primary, primary_header, self.outfile)

            snpeff_anns = annotations[0]
            for snpeff_record in snpeff_anns:
                self.write_record(snpeff_record, self.snpeff_cols, self.snpeff_outfile)

            ma_ann = annotations[1]
            if ma_ann:
                self.write_record(ma_ann, self.ma_cols, self.ma_outfile)

            id_anns = annotations[2]
            for id_ann in id_anns:
                self.write_record(id_ann, self.ids_cols, self.ids_outfile)


if __name__ == "__main__":
    with VcfParser(VCF_FILE, 'parsed', 'snpeff', 'ma', 'ids') as vcf_parser:
        vcf_parser.write()
