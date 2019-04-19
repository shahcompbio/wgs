import os
import pypeliner
from wgs.utils import helpers
import shutil

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
        '--id', sample_id,
        '--WIG', input_wig,
        '--ploidy', "c(2,3)",
        '--normal', "c(0.5,0.6,0.7,0.8,0.9)",
        '--maxCN', 5,
        '--gcWig', gc_wig,
        '--mapWig', map_wig,
        '--normalPanel', normal_panel,
        '--includeHOMD', False,
        '--chrs', "c(1:22, \"X\")",
        '--chrTrain', "c(1:22)",
        '--estimateNormal', True,
        '--estimatePloidy', True,
        '--estimateScPrevalence', True,
        '--scStates', "c(1,3)",
        '--txnE', 0.9999,
        '--txnStrength', 10000,
        '--outDir', plots_dir
    ]

    if centromere:
        cmd.extend(['--centromere', centromere])

    pypeliner.commandline.execute(*cmd)

    segfile = os.path.join(plots_dir, '{}.seg'.format(sample_id))
    shutil.move(segfile, segments)

    paramsfile = os.path.join(plots_dir, '{}.params.txt'.format(sample_id))
    shutil.move(paramsfile, params)

    depthfile = os.path.join(plots_dir, '{}.correctedDepth.txt'.format(sample_id))
    shutil.move(depthfile, depth)

    helpers.make_tarfile(plots_tar, plots_dir)