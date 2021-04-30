import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import alignment


def alignment_metrics_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    assert len(list(inputs.keys())) == 1

    sample_id = list(inputs.keys())[0]
    input_bam = inputs[sample_id]['bam']

    outdir = args['out_dir']

    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    outputs_tdf = os.path.join(outdir, sample_id, '{}.bam.tdf'.format(sample_id))
    metrics_output = os.path.join(outdir, sample_id, '{}_metrics.csv'.format(sample_id))
    metrics_tar = os.path.join(outdir, sample_id, '{}_metrics.tar.gz'.format(sample_id))

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    workflow.subworkflow(
        name="collect_bam_metrics",
        func=alignment.collect_bam_metrics,
        args=(
            mgd.InputFile(input_bam, extensions=['.bai']),
            sample_id,
            args['refdir'],
            mgd.OutputFile(metrics_output),
            mgd.OutputFile(metrics_tar),
            mgd.OutputFile(outputs_tdf),

        )
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            outdir,
            [outputs_tdf,
             metrics_output,
             metrics_tar],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': inputs,
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'alignment'}
        }
    )

    pyp.run(workflow)
