import os
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import ichorcna

def call_cfdna_copynumber(
        samples, config, tumours, normals,
        segments, params, depth, plots_tar
):

    config = config
    segments = dict([(sampid, segments[sampid])
                          for sampid in samples])
    params = dict([(sampid, params[sampid])
                          for sampid in samples])
    depth = dict([(sampid, depth[sampid])
                          for sampid in samples])
    plots_tar = dict([(sampid, plots_tar[sampid])
                          for sampid in samples])


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='ichorcna',
        func=ichorcna.create_ichorcna_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.InputFile('normal_panel', 'sample_id', fnames=normals),
            mgd.OutputFile('segments', 'sample_id', fnames=segments),
            mgd.OutputFile('params', 'sample_id', fnames=params),
            mgd.OutputFile('depth', 'sample_id', fnames=depth),
            config,
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('plots_tar', 'sample_id', fnames=plots_tar),
    ),
    )


    return workflow


def cfdna_cna_calling_workflow(args):

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    config = config['ichorcna']

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}

    cna_outdir = os.path.join(args['out_dir'], 'cfdna_copynumber', '{sample_id}')
    ichor_segments = os.path.join(cna_outdir, 'ichor', 'segments.seg')
    ichor_params = os.path.join(cna_outdir, 'ichor', 'params.txt')
    ichor_corrected_depth = os.path.join(cna_outdir, 'remixt', 'ichor_corrected_depth.txt')
    ichor_plots = os.path.join(cna_outdir, 'remixt', 'plots.tar')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='cfdna_copynumber_calling',
        func=call_cfdna_copynumber,
        args=(
            samples,
            config,
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.panel", 'sample_id', fnames=normals,
                          axes_origin=[]),
            mgd.OutputFile('ichor_segments', 'sample_id', axes_origin=[], template=ichor_segments),
            mgd.OutputFile('ichor_params', 'sample_id', axes_origin=[], template=ichor_params),
            mgd.OutputFile('ichor_corrected_depth', 'sample_id', axes_origin=[], template=ichor_corrected_depth),
            mgd.OutputFile('ichor_plots', 'sample_id', axes_origin=[], template=ichor_plots),
        )
    )

    pyp.run(workflow)
