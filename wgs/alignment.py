import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import alignment


def alignment_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])
    outdir = args['out_dir']

    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs)

    sample_info = helpers.get_sample_info(inputs)

    sample_id = sample_info['SM']

    outputs = os.path.join(outdir, sample_id, '{}.bam'.format(sample_id))
    outputs_tdf = os.path.join(outdir, sample_id, '{}.bam.tdf'.format(sample_id))
    fastqc_tar = os.path.join(outdir, sample_id, '{}_fastqc_metrics.tar.gz'.format(sample_id))
    metrics_output = os.path.join(outdir, sample_id, '{}_metrics.csv'.format(sample_id))
    metrics_tar = os.path.join(outdir, sample_id, '{}_metrics.tar.gz'.format(sample_id))

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    workflow.setobj(
        obj=mgd.OutputChunks('lane_id'),
        value=list(fastqs_r1.keys()),
    )

    workflow.subworkflow(
        name="align_samples",
        func=alignment.align_samples,
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'lane_id', fnames=fastqs_r2),
            mgd.OutputFile(outputs),
            mgd.OutputFile(metrics_output),
            mgd.OutputFile(metrics_tar),
            mgd.OutputFile(fastqc_tar),
            mgd.OutputFile(outputs_tdf),
            sample_info,
            args['refdir'],
            sample_id
        ),
        kwargs={'single_node': args['single_node'],
                'picard_mem': args['picard_mem']}
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            outdir,
            [
                outputs,
                outputs_tdf,
                metrics_output,
                metrics_tar,
            ],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inputs,
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'alignment'}
        }
    )

    pyp.run(workflow)
