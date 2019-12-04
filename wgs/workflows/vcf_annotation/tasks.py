'''
Created on Feb 21, 2018

@author: pwalters
'''
import os

import pandas as pd
import pypeliner
import vcf
import wgs.utils.low_mappability_utils as low_mapp_utils

scripts_directory = os.path.join(
    os.path.realpath(os.path.dirname(__file__)), 'scripts')


def run_snpeff(infile, output, config, docker_image=None):
    '''
    Run snpEff script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    :param config: path to the config YAML file
    '''

    snpeff_config = config['snpeff_params']['snpeff_config']

    cmd = ['snpEff', '-Xmx5G', '-Xms5G', '-XX:ParallelGCThreads=1',
           'GRCh37.75', '-noStats', infile, '>',
           output]

    if snpeff_config:
        cmd.extend(['-c', snpeff_config])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_mutation_assessor(infile, output, config):
    '''
    Run Mutation Assessor script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    mutation_assessor_script = os.path.join(
        scripts_directory, 'annotate_mutation_assessor.py')
    db = config['mutation_assessor_params']['db']

    cmd = ['python', mutation_assessor_script, '--vcf', infile, '--output',
           output, '--db', db]

    pypeliner.commandline.execute(*cmd)


def run_DBSNP(infile, output, config):
    '''
    Run DBSNP script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['dbsnp_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'DBSNP', '--input_type', 'snv', '--flag_with_id']

    pypeliner.commandline.execute(*cmd)


def run_1000gen(infile, output, config):
    '''
    Run 1000Gen script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['thousandgen_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', '1000Gen', '--input_type', 'snv']

    pypeliner.commandline.execute(*cmd)


def run_cosmic(infile, output, config):
    '''
    Run Cosmic script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'add_db_anno.py')
    db = config['cosmic_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'Cosmic', '--input_type', 'snv', '--flag_with_id']

    pypeliner.commandline.execute(*cmd)


def flag_low_mappability(infile, output, blacklist):
    '''
    adds a flag to infile at each row as to
    whether or not they fall within a low-mappability
    region
    :param infile: input vcf
    :param output: output path to vcf with flag
    :param config: config
    '''
    blacklist = pd.read_csv(blacklist)
    vcf_reader = vcf.Reader(open(infile))
    output = vcf.Writer(open(output, "w"), vcf_reader)

    chroms = [record.CHROM for record in vcf_reader]
    positions = [record.POS for record in vcf_reader]

    call_locations = pd.DataFrame({"chromosome": chroms, "positions": positions})

    low_mapp_indexes = low_mapp_utils.is_low_mappability(call_locations, blacklist)

    low_mapp_anno = low_mapp_utils.generate_low_mappability_annotation(low_mapp_indexes,
                                                                       len(vcf_reader))

    for i, record in enumerate(vcf_reader):
        reader.INFO["is_low_mappability"] = low_mapp_anno[i]
        output.write_record(record)


def finalize_vcf(infile, outfile, config):
    # run  cat infile vcf-sort > temp_sorted
    # run  bgzip temp_sorted > outfile (infile + '.gz')
    # run  bcftools index outfile
    # run tabix -f -p vcf outfile
    raise NotImplementedError()
