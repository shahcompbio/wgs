import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from workflows import alignment

from wgs.config import config

def alignment_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])
    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    outputs = os.path.join(outdir, '{sample_id}', '{sample_id}.bam')
    outputs_tdf = os.path.join(outdir, '{sample_id}', '{sample_id}.bam.tdf')
    metrics_output = os.path.join(outdir, '{sample_id}', '{sample_id}_metrics.csv')
    metrics_tar = os.path.join(outdir, '{sample_id}', '{sample_id}_metrics.tar.gz')

    samples = inputs.keys()
    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs, samples, None)

    sample_info = helpers.get_sample_info(inputs)

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane_id'),
        value=fastqs_r1.keys(),
    )

    workflow.subworkflow(
        name="align_samples",
        func=alignment.align_samples,
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r2),
            mgd.Template('output.bam', 'sample_id', template=outputs),
            mgd.Template('metrics.txt', 'sample_id', template=metrics_output),
            mgd.Template('metrics.tar', 'sample_id', template=metrics_tar),
            mgd.Template('output.bam.tdf', 'sample_id', template=outputs_tdf),
            sample_info,
            args['refdir']
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list([outputs, metrics_output, metrics_tar], samples, "sample_id")
    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            outdir,
            outputted_filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inputs,
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'alignment'}
        }
    )

    pyp.run(workflow)
