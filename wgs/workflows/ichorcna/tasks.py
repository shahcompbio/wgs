import os
import sys
import pypeliner
from pypeliner_utils import helpers, pdfutils
import glob

scripts = os.path.join(
    os.path.realpath(os.path.dirname(__file__)),
    'scripts'
)


def hmmcopy_readcounter(input_bam, output_wig, config):
    chromosomes = ['--chromosomes'] + config['hmmcopy_readcounter']['chromosomes']

    cmd = [
              'python',
              os.path.join(scripts, 'read_counter.py'),
              input_bam,
              output_wig,
              '-w',
              str(config['hmmcopy_readcounter']['w']),
              '-m',
              str(config['hmmcopy_readcounter']['m']),
          ] + chromosomes

    helpers.run_cmd(cmd, output=output_wig)


def run_ichorcna(input_wig, output_file, output_obj, config):
    cmd = [
        'Rscript', os.path.join(scripts, 'run_ichorcna.R'),
        input_wig,
        config['calc_corr']['gc'],
        config['calc_corr']['map'],
    ]

    pypeliner.commandline.execute(*cmd)
