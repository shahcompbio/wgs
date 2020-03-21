import scripts.parse_destruct as parse_destruct
import scripts.parse_lumpy as parse_lumpy

from .scripts import consensus


def parse_destruct_task(infile, output, config, chromosomes=None):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''
    if chromosomes:
        config['chromosomes'] = chromosomes

    parse_destruct.parser(infile, output, config, foldback_threshold=config['foldback_threshold'])


def parse_lumpy_task(infile, output, config, chromosomes=None):
    if chromosomes:
        config['chromosomes'] = chromosomes

    vcfdata = parse_lumpy.parse_vcf(infile)
    vcfdata = parse_lumpy.parse_vcf_group(vcfdata)
    vcfdata = parse_lumpy.create_data(vcfdata)
    vcfdata = parse_lumpy.convert_to_df(vcfdata)
    vcfdata = parse_lumpy.filter_calls(vcfdata, config)
    parse_lumpy.write(output, vcfdata)


def consensus_calls(destruct_data, lumpy_data, consensus_calls, config):
    consensus.consensus(
        destruct_data, lumpy_data, consensus_calls, confidence_interval=config['confidence_interval_size']
    )
