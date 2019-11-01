import pypeliner
from wgs.utils import helpers
from workflows import alignment


def alignment_workflow(args):
    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    outputs = helpers.get_values_from_input(inputs, 'bam')
    samples = outputs.keys()
    fastqs_r1, fastqs_r2 = helpers.get_fastqs(inputs, samples, None)

    outdir = args['out_dir']

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
        ),
        kwargs={'single_node': args['single_node']}
    )

    pyp.run(workflow)
