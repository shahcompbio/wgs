import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import remixt
from wgs.workflows import titan


def copynumber_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    global_config = config['globals']
    config = config['copynumber_calling']

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    targets = helpers.get_values_from_input(inputs, 'target_list')
    breakpoints = helpers.get_values_from_input(inputs, 'breakpoints')
    samples = tumours.keys()

    cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')
    remixt_results_filename = os.path.join(cna_outdir, 'remixt', 'results.h5')
    remixt_raw_dir = os.path.join(cna_outdir, 'remixt', 'raw_data')

    titan_raw_dir = os.path.join(cna_outdir, 'titan')
    titan_segments_filename = os.path.join(titan_raw_dir, 'segments.h5')
    titan_markers_filename = os.path.join(titan_raw_dir, 'markers.h5')
    titan_params_filename = os.path.join(titan_raw_dir, 'params.h5')

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs'])
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='titan',
        func=titan.create_titan_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("target_list", 'sample_id', fnames=targets,
                          axes_origin=[]),
            mgd.Template(titan_raw_dir, 'sample_id'),
            mgd.OutputFile('titan_segments_filename', 'sample_id',
                           axes_origin=[], template=titan_segments_filename),
            mgd.OutputFile('titan_params_filename', 'sample_id',
                           axes_origin=[], template=titan_params_filename),
            mgd.OutputFile('titan_markers_filename', 'sample_id',
                           axes_origin=[], template=titan_markers_filename),
            global_config,
            config,
            config['titan_intervals'],
            mgd.InputInstance('sample_id'),
        ),
        kwargs={'single_node': args['single_node']}
    )

    workflow.subworkflow(
        name='remixt',
        func=remixt.create_remixt_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_bam', 'sample_id',
                          fnames=tumours, extensions=['.bai']),
            mgd.InputFile('normal_bam', 'sample_id',
                          fnames=normals, extensions=['.bai']),
            mgd.InputFile('destruct_breakpoints', 'sample_id',
                          axes_origin=[], fnames=breakpoints),
            mgd.InputInstance('sample_id'),
            config['remixt_refdata'],
            mgd.OutputFile('remixt_results_filename', 'sample_id',
                           axes_origin=[], template=remixt_results_filename),
            mgd.Template(remixt_raw_dir, 'sample_id'),
            config['min_num_reads'],
            global_config
        ),
        kwargs={'single_node': args['single_node'],
                'docker_containers': config['docker']}
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            [titan_segments_filename, titan_params_filename,
             titan_markers_filename, remixt_results_filename],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'realignment'}
        }
    )

    pyp.run(workflow)