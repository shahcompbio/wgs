import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.utils import helpers

from wgs.workflows import realignment


def realign_bams(samples, inputs, outputs, outdir, config, single_node=False):
    outputs = dict([(sampid, outputs[sampid])
                    for sampid in samples])
    inputs = dict([(sampid, inputs[sampid])
                    for sampid in samples])


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='realign_bam_file',
        func=realignment.realign_bam_file,
        axes=('sample_id',),
        args=(
            mgd.InputFile("input.bam", "sample_id", fnames=inputs),
            mgd.OutputFile("output.bam", "sample_id", fnames=outputs),
            outdir,
            config
        ),
        kwargs={'single_node': single_node}
    )

    return workflow


def realign_bam_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    config = config['alignment']

    yamldata = yaml.safe_load(open(args['input_yaml']))

    samples = yamldata.keys()

    input_bams = {sample: yamldata[sample]['input'] for sample in samples}
    output_bams = {sample: yamldata[sample]['output'] for sample in samples}

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
            mgd.OutputFile("realigned.bam", 'sample_id', fnames=output_bams,
                           extensions=['.bai'], axes_origin=[]),
            args['out_dir'],
            config
        ),
        kwargs={'single_node': args['single_node']}
    )

    pyp.run(workflow)
