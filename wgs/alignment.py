import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from workflows import alignment


def paired_alignment(
        config, tumours, normals, samples, tumour_fastqs_r1,
        tumour_fastqs_r2, normal_fastqs_r1, normal_fastqs_r2,
        outdir,
        single_node=False
):
    tumours = dict([(sampid, tumours[sampid])
                    for sampid in samples])
    normals = dict([(sampid, normals[sampid])
                    for sampid in samples])

    workflow = pypeliner.workflow.Workflow()

    global_config = config['globals']
    config = config['alignment']

    tumour_outdir = os.path.join(outdir, 'tumour')
    normal_outdir = os.path.join(outdir, 'normal')

    workflow.setobj(
        obj=mgd.OutputChunks('tum_sample_id', 'tum_lane'),
        value=tumour_fastqs_r1.keys(),
    )

    workflow.setobj(
        obj=mgd.OutputChunks('norm_sample_id', 'norm_lane'),
        value=normal_fastqs_r1.keys(),
    )

    workflow.subworkflow(
        name='align_tumours',
        func=alignment.align_samples,
        args=(
            config,
            global_config,
            tumour_fastqs_r1,
            tumour_fastqs_r2,
            mgd.OutputFile('tumour.bam', 'tum_sample_id', fnames=tumours, extensions=['.bai'], axes_origin=[]),
            tumour_outdir
        )
    )

    workflow.subworkflow(
        name='align_normals',
        func=alignment.align_samples,
        args=(
            config,
            global_config,
            normal_fastqs_r1,
            normal_fastqs_r2,
            mgd.OutputFile('normal.bam', 'tum_sample_id', fnames=normals, extensions=['.bai'], axes_origin=[]),
            normal_outdir
        )
    )

    return workflow


def alignment_workflow(args):
    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    outputs = helpers.get_values_from_input(inputs, 'bam')
    samples = outputs.keys()
    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs, samples, None)

    outdir = args['out_dir']

    config_globals = config['globals']
    config = config['alignment']

    ctx = {'docker_image': config['docker']['wgs']}

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

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
        ),
        kwargs={'single_node': args['single_node']}
    )

    pyp.run(workflow)
