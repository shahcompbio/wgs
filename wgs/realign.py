import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import realignment

def realign_bam_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    input_bam = args['input_bam']
    sample_id = args['sample_id']

    output_bam = os.path.join(outdir, sample_id, '{}.bam'.format(sample_id))
    metrics = os.path.join(outdir, sample_id, '{}_metrics.csv'.format(sample_id))
    metrics_tar = os.path.join(outdir, sample_id, '{}.tar'.format(sample_id))

    workflow.subworkflow(
        name="realign",
        func=realignment.realign_bam_files,
        ctx=helpers.get_default_ctx(),
        args=(
            mgd.InputFile(input_bam, extensions=['.bai']),
            mgd.OutputFile(output_bam, extensions=['.bai']),
            mgd.OutputFile(metrics),
            mgd.OutputFile(metrics_tar),
            args['refdir'],
        ),
        kwargs={
            'single_node': args['single_node'],
            'ignore_bamtofastq_exception': args['ignore_bamtofastq_exception'],
            'picard_mem': args['picard_mem']
        }
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            [output_bam, metrics, metrics_tar],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'realignment'}
        }
    )

    pyp.run(workflow)
