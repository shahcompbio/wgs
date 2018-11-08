'''
Created on Feb 21, 2018

@author: pwalters
'''
import os
import pypeliner

scripts_directory = os.path.join(
    os.path.realpath(os.path.dirname(__file__)), 'scripts')

def run_snpeff(infile, output, config):
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

    pypeliner.commandline.execute(*cmd)


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

    script = os.path.join(scripts_directory, 'flagPos.py')
    db = config['dbsnp_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'DBSNP', '--input_type', 'snv', '--flag_with_id', 'True']

    pypeliner.commandline.execute(*cmd)


def run_1000gen(infile, output, config):
    '''
    Run 1000Gen script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'flagPos.py')
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

    script = os.path.join(scripts_directory, 'flagPos.py')
    db = config['cosmic_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'Cosmic', '--input_type', 'snv', '--flag_with_id',
           'True']

    pypeliner.commandline.execute(*cmd)

def finalize_vcf(infile, outfile, config):

     raise NotImplementedError()

    # run  cat infile vcf-sort > temp_sorted
    # run  bgzip temp_sorted > outfile (infile + '.gz')
    # run  bcftools index outfile
    # run tabix -f -p vcf outfile
