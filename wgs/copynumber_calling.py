import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import remixt
from wgs.workflows import titan


def copynumber_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    run_titan = args['titan']
    run_remixt = args['remixt']

    if not run_titan and not run_remixt:
        run_titan = True
        run_remixt = True

    inputs = helpers.load_yaml(args['input_yaml'])

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    targets = helpers.get_values_from_input(inputs, 'target_list')
    breakpoints = helpers.get_values_from_input(inputs, 'breakpoints')
    samples = list(tumours.keys())

    cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')

    titan_raw_dir = os.path.join(cna_outdir, 'titan')

    titan_outfile = os.path.join(titan_raw_dir, '{sample_id}_titan_markers.csv.gz')
    titan_params = os.path.join(titan_raw_dir, '{sample_id}_titan_params.csv.gz')
    titan_segs = os.path.join(titan_raw_dir, '{sample_id}_titan_segs.csv.gz')
    titan_igv_segs = os.path.join(titan_raw_dir, '{sample_id}_titan_igv_segs.seg')
    titan_parsed = os.path.join(titan_raw_dir, '{sample_id}_titan_parsed.csv.gz')
    titan_plots = os.path.join(titan_raw_dir, '{sample_id}_titan_plots.pdf')
    titan_tar_outputs = os.path.join(titan_raw_dir, '{sample_id}_data_all_parameters.tar.gz')
    museq_vcf = os.path.join(titan_raw_dir, '{sample_id}_museq.vcf')

    remixt_outdir = os.path.join(args['out_dir'], 'remixt', '{sample_id}')
    remixt_outfile = os.path.join(remixt_outdir, '{sample_id}_remixt.h5')
    remixt_raw_dir = os.path.join(remixt_outdir, '{sample_id}_raw_dir')

    remixt_brk_cn_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_brk_cn.csv.gz')
    remixt_cn_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_cn.csv.gz')
    remixt_minor_modes_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_minor_modes.csv.gz')
    remixt_mix_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_mix.csv.gz')
    remixt_read_depth_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_read_depth.csv.gz')
    remixt_stats_csv = os.path.join(remixt_outdir, '{sample_id}_remixt_stats.csv.gz')

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    if run_remixt:
        workflow.subworkflow(
            name='remixt',
            func=remixt.create_remixt_workflow,
            axes=('sample_id',),
            args=(
                mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                              extensions=['.bai']),
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai']),
                mgd.InputFile("breakpoints", 'sample_id', fnames=breakpoints),
                mgd.InputInstance('sample_id'),
                mgd.OutputFile('remixt.h5', 'sample_id', template=remixt_outfile),
                mgd.OutputFile('remixt_brk_cn.csv', 'sample_id', template=remixt_brk_cn_csv),
                mgd.OutputFile('remixt_cn.csv', 'sample_id', template=remixt_cn_csv),
                mgd.OutputFile('remixt_minor_modes.csv', 'sample_id', template=remixt_minor_modes_csv),
                mgd.OutputFile('remixt_mix.csv', 'sample_id', template=remixt_mix_csv),
                mgd.OutputFile('remixt_read_depth.csv', 'sample_id', template=remixt_read_depth_csv),
                mgd.OutputFile('remixt_stats.csv', 'sample_id', template=remixt_stats_csv),
                refdir_paths['refdata_remixt'],
                mgd.Template('rawdir', 'sample_id', template=remixt_raw_dir),
                refdir_paths['reference'],
                chromosomes
            ),
            kwargs={'single_node': args['single_node']}
        )

    if run_titan:
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
                mgd.OutputFile('museq.vcf', 'sample_id', template=museq_vcf),
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

    filenames = []

    if run_remixt:
        filenames += [
            remixt_outfile,
            remixt_raw_dir,
            remixt_brk_cn_csv,
            remixt_cn_csv,
            remixt_minor_modes_csv,
            remixt_mix_csv,
            remixt_read_depth_csv,
            remixt_stats_csv
        ]
    if run_titan:
        filenames += [
            titan_outfile,
            titan_params,
            titan_segs,
            titan_igv_segs,
            titan_parsed,
            titan_plots,
            titan_tar_outputs,
            museq_vcf,
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
            'metadata': {'type': 'copynumber_calling'}
        }
    )

    pyp.run(workflow)
