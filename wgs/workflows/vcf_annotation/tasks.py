"""
Created on Feb 21, 2018

@author: pwalters
"""
import os

import pypeliner
from wgs.utils import helpers
from wgs.workflows.vcf_annotation.scripts import flag_mappability

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')


def split_by_chrom(inputfile, outputfile, chrom):
    with helpers.GetFileHandle(inputfile) as infile, helpers.GetFileHandle(outputfile, 'wt') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                if line.split('\t')[0] == chrom:
                    outfile.write(line)


def run_snpeff(infile, output, config, docker_image=None):
    """
    Run snpEff script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    :param config: path to the config YAML file
    """
    snpeff_config = config['snpeff_params']['snpeff_config']
    cmd = [
        'snpEff', '-Xmx5G', '-Xms5G', '-XX:ParallelGCThreads=1',
        'GRCh37.75', '-noStats', infile, '>',
        output]
    if snpeff_config:
        cmd.extend(['-c', snpeff_config])
    pypeliner.commandline.execute(docker_image=docker_image, *cmd)


def run_mutation_assessor(infile, output, config, chrom):
    """
    Run Mutation Assessor script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    """
    mutation_assessor_script = os.path.join(scripts_directory, 'annotate_mutation_assessor.py')
    db = config['mutation_assessor_params']['db']
    cmd = [
        'python', mutation_assessor_script, '--vcf', infile, '--output',
        output, '--db', db, '--chrom', chrom]
    pypeliner.commandline.execute(*cmd)


def run_DBSNP(infile, output, config, chrom, input_type='snv'):
    """
    Run DBSNP script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    """
    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['dbsnp_params']['db']
    cmd = [
        'python', script, '--infile', infile, '--db', db,
        '--label', 'DBSNP', '--input_type', input_type, '--flag_with_id',
        '--out', output, '--chrom', chrom]

    pypeliner.commandline.execute(*cmd)


def run_1000gen(infile, output, config, chrom, input_type='snv'):
    """
    Run 1000Gen script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    """
    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['thousandgen_params']['db']
    cmd = [
        'python', script, '--infile', infile, '--db', db,
        '--label', '1000Gen', '--input_type', input_type,
        '--out', output, '--chrom', chrom]

    pypeliner.commandline.execute(*cmd)


def run_cosmic(infile, output, config, chrom, input_type='snv'):
    """
    Run Cosmic script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    """
    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['cosmic_params']['db']
    cmd = [
        'python', script, '--infile', infile, '--db', db,
        '--label', 'Cosmic', '--input_type', input_type, '--flag_with_id',
        '--out', output, '--chrom', chrom]

    pypeliner.commandline.execute(*cmd)


def flag_low_mappability(infile, output, blacklist, chrom):
    """
    adds a flag to infile at each row as to
    whether or not they fall within a low-mappability
    region
    :param infile: input vcf
    :param output: output path to vcf with flag
    :param config: config
    """
    flag_mappability.main(infile, output, blacklist, chrom)
