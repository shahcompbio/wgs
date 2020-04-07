import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import realignment


def realign_bams(samples, inputs, outputs, metrics, metrics_tar, refdir, single_node=False):
    outputs = dict([(sampid, outputs[sampid])
                    for sampid in samples])
    inputs = dict([(sampid, inputs[sampid])
                   for sampid in samples])

    metrics = dict([(sampid, metrics[sampid])
                    for sampid in samples])
    metrics_tar = dict([(sampid, metrics_tar[sampid])
                        for sampid in samples])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='realign_bam_file',
        func=realignment.realign_bam_files,
        args=(
            mgd.InputFile("input.bam", "sample_id", axes_origin=[], fnames=inputs),
            mgd.OutputFile("output.bam", "sample_id", axes_origin=[], fnames=outputs),
            mgd.OutputFile("output.txt", "sample_id", axes_origin=[], fnames=metrics),
            mgd.OutputFile("output.tar", "sample_id", axes_origin=[], fnames=metrics_tar),
            refdir,
            samples
        ),
        kwargs={'single_node': single_node}
    )

    return workflow


def realign_bam_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    yamldata = yaml.safe_load(open(args['input_yaml']))

    samples = list(yamldata.keys())

    input_bams = {sample: yamldata[sample]['input'] for sample in samples}

    output_bams = os.path.join(outdir, '{sample_id}', '{sample_id}.bam')
    metrics = os.path.join(outdir, '{sample_id}', '{sample_id}.txt')
    metrics_tar = os.path.join(outdir, '{sample_id}', '{sample_id}.tar')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="realign",
        func=realign_bams,
        ctx=helpers.get_default_ctx(),
        args=(
            samples,
            mgd.InputFile("input.bam", 'sample_id', fnames=input_bams,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile("realigned.bam", 'sample_id', template=output_bams,
                           extensions=['.bai', '.tdf'], axes_origin=[]),
            mgd.OutputFile("realigned.txt", 'sample_id', template=metrics,
                           extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile("realigned.tar", 'sample_id', template=metrics_tar,
                           extensions=['.bai'], axes_origin=[]),
            args['refdir'],
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list([output_bams, metrics, metrics_tar], samples, 'sample_id')

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
