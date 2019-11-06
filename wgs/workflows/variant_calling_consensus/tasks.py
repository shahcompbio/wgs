from wgs.utils import csvutils

from .scripts import vcfparser


def parse_vcf(infile, primary_table, snpeff_table, ma_table, id_table, parse_config, ):
    with vcfparser.VcfParser(infile, primary_table, snpeff_table, ma_table, id_table) as vcf_parser:
        vcf_parser.write()


def merge_overlap(infiles, outfile, on=('chrom', 'pos', 'ref', 'alt')):
    csvutils.merge_csv(
        infiles, outfile, 'inner', on,
    )
