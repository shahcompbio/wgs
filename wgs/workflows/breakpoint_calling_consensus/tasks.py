from .scripts import consensus
from .scripts import parse_destruct
from .scripts import parse_lumpy


def parse_destruct(infile, output, config):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''

    parse_destruct.parser(infile, output, foldback_threshold=config['foldback_threshold'])


def parse_lumpy(infile, output, config):
    vcfdata = parse_lumpy.parse_vcf(infile)
    parse_lumpy.write(output, vcfdata)


def consensus_calls(destruct_data, lumpy_data, consensus_calls, config):
    consensus.consensus(
        destruct_data, lumpy_data, consensus_calls, confidence_interval=config['confidence_interval']
    )
