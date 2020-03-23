import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import hmmcopy
from wgs.workflows import titan


from wgs.config import config

def copynumber_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    inputs = helpers.load_yaml(args['input_yaml'])

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    targets = helpers.get_values_from_input(inputs, 'target_list')
    samples = tumours.keys()

    cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')

    titan_raw_dir = os.path.join(cna_outdir, 'titan')

    titan_outfile = os.path.join(titan_raw_dir, '{sample_id}_titan_markers.csv.gz')
    titan_params = os.path.join(titan_raw_dir, '{sample_id}_titan_params.csv.gz')
    titan_segs = os.path.join(titan_raw_dir, '{sample_id}_titan_segs.csv.gz')
    titan_igv_segs = os.path.join(titan_raw_dir, '{sample_id}_titan_igv_segs.seg')
    titan_parsed = os.path.join(titan_raw_dir, '{sample_id}_titan_parsed.csv.gz')
    titan_plots = os.path.join(titan_raw_dir, '{sample_id}_titan_plots.pdf')
    titan_tar_outputs = os.path.join(titan_raw_dir, '{sample_id}_data_all_parameters.tar.gz')

    hmmcopy_normal_raw_dir = os.path.join(cna_outdir, 'hmmcopy_normal')
    normal_bias_pdf = os.path.join(hmmcopy_normal_raw_dir, 'plots', '{sample_id}_bias.pdf')
    normal_correction_pdf = os.path.join(hmmcopy_normal_raw_dir, 'plots', '{sample_id}_correction.pdf')
    normal_hmmcopy_pdf = os.path.join(hmmcopy_normal_raw_dir, 'plots', '{sample_id}_hmmcopy.pdf')
    normal_correction_table = os.path.join(hmmcopy_normal_raw_dir, '{sample_id}_correctreads_with_state.txt')
    normal_pygenes = os.path.join(hmmcopy_normal_raw_dir, '{sample_id}_hmmcopy.seg.pygenes')

    hmmcopy_tumour_raw_dir = os.path.join(cna_outdir, 'hmmcopy_tumour')
    tumour_bias_pdf = os.path.join(hmmcopy_tumour_raw_dir, 'plots', '{sample_id}_bias.pdf')
    tumour_correction_pdf = os.path.join(hmmcopy_tumour_raw_dir, 'plots', '{sample_id}_correction.pdf')
    tumour_hmmcopy_pdf = os.path.join(hmmcopy_tumour_raw_dir, 'plots', '{sample_id}_hmmcopy.pdf')
    tumour_correction_table = os.path.join(hmmcopy_tumour_raw_dir, '{sample_id}_correctreads_with_state.txt')
    tumour_pygenes = os.path.join(hmmcopy_tumour_raw_dir, '{sample_id}_hmmcopy.seg.pygenes')

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']


    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
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
                          extensions=['.bai']),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai']),
            mgd.InputFile("target_list", 'sample_id', fnames=targets),
            mgd.OutputFile('outfile', 'sample_id', template=titan_outfile),
            mgd.OutputFile('params', 'sample_id', template=titan_params),
            mgd.OutputFile('segs', 'sample_id', template=titan_segs),
            mgd.OutputFile('igv_segs', 'sample_id', template=titan_igv_segs),
            mgd.OutputFile('parsed', 'sample_id', template=titan_parsed),
            mgd.OutputFile('plots', 'sample_id', template=titan_plots),
            mgd.OutputFile('tar_outputs', 'sample_id', template=titan_tar_outputs),
            mgd.Template(titan_raw_dir, 'sample_id'),
            mgd.InputInstance('sample_id'),
            refdir_paths['reference'],
            chromosomes,
            refdir_paths['het_positions_titan'],
            refdir_paths['map_wig'],
            refdir_paths['gc_wig'],
            refdir_paths['gtf'],
        ),
        kwargs={'single_node': args['single_node']}
    )

    workflow.subworkflow(
        name='hmmcopy_normal',
        func=hmmcopy.create_hmmcopy_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai']),
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('normal_bias', 'sample_id', template=normal_bias_pdf),
            mgd.OutputFile('normal_correction', 'sample_id', template=normal_correction_pdf),
            mgd.OutputFile('normal_hmmcopy', 'sample_id', template=normal_hmmcopy_pdf),
            mgd.OutputFile('normal_correction_table', 'sample_id', template=normal_correction_table),
            mgd.OutputFile('normal_pygenes', 'sample_id', template=normal_pygenes),
            chromosomes,
            refdir_paths['map_wig'],
            refdir_paths['gc_wig'],
            refdir_paths['gtf']
        ),
    )

    workflow.subworkflow(
        name='hmmcopy_tumour',
        func=hmmcopy.create_hmmcopy_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai']),
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('tumour_bias', 'sample_id', template=tumour_bias_pdf),
            mgd.OutputFile('tumour_correction', 'sample_id', template=tumour_correction_pdf),
            mgd.OutputFile('tumour_hmmcopy', 'sample_id', template=tumour_hmmcopy_pdf),
            mgd.OutputFile('tumour_correction_table', 'sample_id', template=tumour_correction_table),
            mgd.OutputFile('tumour_pygenes', 'sample_id', template=tumour_pygenes),
            chromosomes,
            refdir_paths['map_wig'],
            refdir_paths['gc_wig'],
            refdir_paths['gtf']
        ),
    )

    filenames = [
        titan_outfile,
        titan_params,
        titan_segs,
        titan_igv_segs,
        titan_parsed,
        titan_plots,
        titan_tar_outputs,
    ]

    outputted_filenames = helpers.expand_list(filenames, samples, "sample_id")

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            outputted_filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'realignment'}
        }
    )

    pyp.run(workflow)
