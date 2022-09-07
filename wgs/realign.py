import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.utils import helpers
from wgs.workflows import realignment


def realign_bams(
        input, output, metrics,
        metrics_tar, refdir, ignore_bamtofastq_exception,
        single_node=False, picard_mem=8
):

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='realign_bam_file',
        func=realignment.realign_bam_files,
        args=(
            mgd.InputFile(input),
            mgd.OutputFile(output),
            mgd.OutputFile(metrics),
            mgd.OutputFile(metrics_tar),
            refdir,
        ),
        kwargs={
            'single_node': single_node,
            'ignore_bamtofastq_exception': ignore_bamtofastq_exception,
            'picard_mem': picard_mem
        }
    )

    return workflow


def realign_bam_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    yamldata = yaml.safe_load(open(args['input_yaml']))

    input_bam = yamldata['input']

    output_bams = args['output_prefix'] + '_realigned.bam'
    metrics = args['output_prefix'] + '_realigned_metrics.csv'
    metrics_tar = args['output_prefix'] + '_realigned.tar'

    workflow.subworkflow(
        name="realign",
        func=realign_bams,
        ctx=helpers.get_default_ctx(),
        args=(
            mgd.InputFile(input_bam, extensions=['.bai']),
            mgd.OutputFile(output_bams, extensions=['.bai', '.tdf']),
            mgd.OutputFile(metrics, extensions=['.bai']),
            mgd.OutputFile(metrics_tar, extensions=['.bai']),
            args['refdir'],
            args['ignore_bamtofastq_exception']
        ),
        kwargs={
            'single_node': args['single_node'],
            'picard_mem': args['picard_mem']
        }
    )

    outputted_filenames = [output_bams, metrics, metrics_tar]

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            outputted_filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'realignment'}
        }
    )

    pyp.run(workflow)
