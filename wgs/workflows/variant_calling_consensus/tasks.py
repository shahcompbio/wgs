from wgs.utils import csvutils
from wgs.utils import helpers
import os

from .scripts import vcfparser


def parse_vcf(
        infile, primary_table, snpeff_table,
        ma_table, id_table, parse_config, chromosomes, tempdir
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

    helpers.makedirs(tempdir)

    primary_temp = os.path.join(tempdir, 'primary.csv')
    snpeff_temp = os.path.join(tempdir, 'snpeff.csv')
    ma_temp = os.path.join(tempdir, 'ma.csv')
    ids_temp = os.path.join(tempdir, 'ids.csv')


    filter_out = []
    if 'filter_low_mappability' in parse_config and parse_config['filter_low_mappability']:
        filter_out.append(('LOW_MAPPABILITY', 'eq', True))

    if chromosomes:
        filter_out.append(('CHROM', 'notin', chromosomes))

    if 'pr_threshold' in parse_config and parse_config['pr_threshold']:
        filter_out.append(('PR', 'lt', parse_config['pr_threshold']))

    with vcfparser.VcfParser(
            infile, primary_temp, snpeff_temp,
            ma_temp, ids_temp, filter_out) as vcf_parser:
        vcf_parser.write()

    csvutils.finalize_csv(primary_temp, primary_table)
    csvutils.finalize_csv(snpeff_temp, snpeff_table)
    csvutils.finalize_csv(ma_temp, ma_table)
    csvutils.finalize_csv(ids_temp, id_table)


def merge_overlap(infiles, outfile, on=('chrom', 'pos', 'ref', 'alt')):
    csvutils.merge_csv(
        infiles, outfile, 'inner', on,
    )
