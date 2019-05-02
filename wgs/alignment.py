import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers

from workflows import alignment


def paired_alignment(
        config, tumours, normals, samples, tumour_fastqs_r1,
        tumour_fastqs_r2, normal_fastqs_r1, normal_fastqs_r2,
        outdir_template_normal, outdir_template_tumour,
        single_node=False
):
    tumours = dict([(sampid, tumours[sampid])
                    for sampid in samples])
    normals = dict([(sampid, normals[sampid])
                    for sampid in samples])
    tumours_index = dict([(sampid, tumours[sampid] + '.bai')
                          for sampid in samples])
    normals_index = dict([(sampid, normals[sampid] + '.bai')
                          for sampid in samples])

    workflow = pypeliner.workflow.Workflow()

    global_config = config['globals']
    config = config['alignment']

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
        func=alignment.align_sample,
        axes=('tum_sample_id', 'tum_lane'),
        args=(
            config,
            mgd.InputFile('input.r1.fastq.gz', 'tum_sample_id', 'tum_lane', fnames=tumour_fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'tum_sample_id', 'tum_lane', fnames=tumour_fastqs_r2),
            mgd.TempOutputFile('tumour.bam', 'tum_sample_id', 'tum_lane'),
            mgd.Template(outdir_template_tumour, 'tum_sample_id', 'tum_lane'),
            [mgd.InputInstance('tum_sample_id'),
             mgd.InputInstance('tum_lane')]
        ),
        kwargs={'single_node': single_node}
    )

    workflow.transform(
        name='merge_tumour_lanes',
        ctx={'mem': global_config['memory']['med'], 'ncpus': 1},
        func="wgs.workflows.alignment.tasks.merge_bams",
        axes=('tum_sample_id',),
        args=(
            mgd.TempInputFile('tumour.bam', 'tum_sample_id', 'tum_lane'),
            mgd.OutputFile('output.bam', 'tum_sample_id', fnames=tumours),
            None
        )
    )

    workflow.subworkflow(
        name='align_normals',
        func=alignment.align_sample,
        axes=('norm_sample_id', 'norm_lane'),
        args=(
            config,
            mgd.InputFile('input.r1.fastq.gz', 'norm_sample_id', 'norm_lane', fnames=normal_fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'norm_sample_id', 'norm_lane', fnames=normal_fastqs_r2),
            mgd.TempOutputFile('normal.bam', 'norm_sample_id', 'norm_lane'),
            mgd.Template(outdir_template_normal, 'norm_sample_id', 'norm_lane'),
            [mgd.InputInstance('norm_sample_id'),
             mgd.InputInstance('norm_lane')]
        ),
        kwargs={'single_node': single_node}
    )

    workflow.transform(
        name='merge_normal_lanes',
        ctx={'mem': global_config['memory']['med'], 'ncpus': 1},
        func="wgs.workflows.alignment.tasks.merge_bams",
        axes=('norm_sample_id',),
        args=(
            mgd.TempInputFile('normal.bam', 'norm_sample_id', 'norm_lane'),
            mgd.OutputFile('output.bam', 'norm_sample_id', fnames=normals),
            None
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
