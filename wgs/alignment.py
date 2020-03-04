import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from workflows import alignment


def alignment_workflow(args):
    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])
    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    outputs = helpers.get_values_from_input(inputs, 'bam')
    samples = outputs.keys()
    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs, samples, None)

    sample_info = helpers.get_sample_info(inputs)

    config_globals = config['globals']
    config = config['alignment']

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs']))

    workflow.subworkflow(
        name="align_samples",
        func=alignment.align_samples,
        args=(
            config,
            config_globals,
            fastqs_r1,
            fastqs_r2,
            outputs,
            outdir,
            sample_info,
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list(outputs, samples, 'sample_id')

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
