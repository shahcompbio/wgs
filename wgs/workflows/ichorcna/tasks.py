import os
import sys
import pypeliner
from wgs.utils import helpers
import shutil


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


def run_ichorcna(input_wig, normal_panel, segments,
                 params, depth, centromere,
                 gc_wig, map_wig, sample_id, plots_dir,
                 plots_tar,
                 txnE=None, chromosomes=None):
    cmd = [
        'Rscript', os.path.join(scripts, 'run_ichorcna.R'),
        input_wig,
        normal_panel,
        centromere,
        gc_wig,
        map_wig,
        sample_id,
        txnE,
        '--outDir',
        plots_dir
    ]

    pypeliner.commandline.execute(*cmd)

    segfile = os.path.join(plots_dir, '{}.seg'.format(sample_id))
    shutil.move(segfile, segments)

    paramsfile = os.path.join(plots_dir, '{}.params.txt'.format(sample_id))
    shutil.move(paramsfile, params)

    depthfile = os.path.join(plots_dir, '{}.correctedDepth.txt'.format(sample_id))
    shutil.move(depthfile, depth)

    helpers.make_tarfile(plots_tar, plots_dir)