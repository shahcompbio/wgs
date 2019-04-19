import os
import pypeliner
from wgs.utils import helpers
import shutil
from scripts import ReadCounter


def hmmcopy_readcounter(input_bam, output_wig, config):

    rc = ReadCounter(
        input_bam, output_wig, config['hmmcopy_readcounter']['w'],
        config['chromosomes'], config['hmmcopy_readcounter']['m'],
    )
    rc.main()

def run_ichorcna(
        input_wig, normal_panel, segments, params, depth, centromere, gc_wig,
        map_wig, sample_id, plots_dir, plots_tar, **kwargs
):

    cmd = [
        'runIchorCNA.R',
        '--id', sample_id,
        '--WIG', input_wig,
        '--gcWig', gc_wig,
        '--mapWig', map_wig,
        '--normalPanel', normal_panel,
        '--outDir', plots_dir
    ]

    if centromere:
        cmd.extend(['--centromere', centromere])

    for flag, value in kwargs.items():
        if isinstance(value, list):
            value = 'c(' + ','.join(map(str,value)) + ')'
        elif isinstance(value, bool):
            value = 'True' if value else 'False'

        cmd.extend(['--{}'.format(flag), value])

    pypeliner.commandline.execute(*cmd)

    segfile = os.path.join(plots_dir, '{}.seg'.format(sample_id))
    shutil.move(segfile, segments)

    paramsfile = os.path.join(plots_dir, '{}.params.txt'.format(sample_id))
    shutil.move(paramsfile, params)

    depthfile = os.path.join(plots_dir, '{}.correctedDepth.txt'.format(sample_id))
    shutil.move(depthfile, depth)

    helpers.make_tarfile(plots_tar, plots_dir)