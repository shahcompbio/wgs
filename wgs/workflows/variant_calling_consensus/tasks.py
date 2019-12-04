from wgs.utils import csvutils

from .scripts import vcfparser


def parse_vcf(
        infile, primary_table, snpeff_table,
        ma_table, id_table, parse_config
):
    '''
    parses a vcf containing variant calls
    to a CSV.

    :param infile: vcf containing calls and annotations
    :param primary_table:csv output filepath containing base calls
    :param snpeff_table: csv output filepath containing snpeff annotations
    :param ma_table: csv output filepath containing ma annotations
    :param id_table: csv output filepath containg id annotations
    :param parse_config: config?? currently unused
    :param parse_low_mappability: boolean; whether or not to filter low-mappability calls
    ##assuming there will by a path to a blacklisted calls table in config
    '''

    parse_low_mappability = parse_config['parse_low_mappability']

    with vcfparser.VcfParser(infile, primary_table, snpeff_table,
                             ma_table, id_table, parse_low_mappability) as vcf_parser:
        vcf_parser.write()

    # add a yaml ext to parsed files
    csvutils.prep_csv_files(primary_table, primary_table)
    csvutils.prep_csv_files(snpeff_table, snpeff_table)
    csvutils.prep_csv_files(ma_table, ma_table)
    csvutils.prep_csv_files(id_table, id_table)


def merge_overlap(infiles, outfile, on=('chrom', 'pos', 'ref', 'alt')):
    csvutils.merge_csv(
        infiles, outfile, 'inner', on,
    )
