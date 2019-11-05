from wgs.utils import csvutils

from .scripts import vcfparser


def parse_museq(infile, primary_table, snpeff_table, parse_config):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''
    vcfdata = vcfparser.parse_vcf(infile)
    vcfparser.write({'primary': primary_table, 'snpeff': snpeff_table}, vcfdata)


def parse_strelka(infile, primary_table, snpeff_table, parser_config):
    vcfdata = vcfparser.parse_vcf(infile)
    vcfparser.write({'primary': primary_table, 'snpeff': snpeff_table}, vcfdata)


def merge_overlap(infiles, outfile):
    csvutils.merge_csv(
        infiles, outfile, 'inner', ['chrom', 'pos', 'ref', 'alt'],
    )
