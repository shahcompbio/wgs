import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import alignment


def alignment_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    outputs = args['output_prefix'] + 'aligned.bam'
    outputs_tdf = args['output_prefix'] + 'aligned.bam.tdf'
    metrics_output = args['output_prefix'] + 'aligned_metrics.csv'
    metrics_tar = args['output_prefix'] + 'aligned_metrics.tar.gz'

    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs)

    sample_info = inputs['readgroup_info']

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

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
            mgd.Template(outputs),
            mgd.Template(metrics_output),
            mgd.Template(metrics_tar),
            mgd.Template(outputs_tdf),
            sample_info,
            args['refdir']
        ),
        kwargs={'single_node': args['single_node'],
                'picard_mem': args['picard_mem']}
    )

    outputted_filenames = [
        outputs,
        outputs_tdf,
        metrics_output,
        metrics_tar,

    ]

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
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
