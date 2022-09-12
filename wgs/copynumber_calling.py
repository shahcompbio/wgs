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

    sex = args['sex']

    if not run_titan and not run_remixt:
        run_titan = True
        run_remixt = True

    inputs = helpers.load_yaml(args['input_yaml'])

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    normal = inputs['normal']
    tumour = inputs['tumour']
    targets = inputs['target_list']
    sample_id = inputs['sample_id']

    titan_outfile = args['output_prefix'] + '_titan_markers.csv.gz'
    titan_params = args['output_prefix'] + '_titan_params.csv.gz'
    titan_segs = args['output_prefix'] + '_titan_segs.csv.gz'
    titan_igv_segs = args['output_prefix'] + '_titan_igv_segs.seg'
    titan_parsed = args['output_prefix'] + '_titan_parsed.csv.gz'
    titan_plots = args['output_prefix'] + '_titan_plots.pdf'
    titan_tar_outputs = args['output_prefix'] + '_data_all_parameters.tar.gz'
    museq_vcf = args['output_prefix'] + '_museq.vcf'

    remixt_outfile = args['output_prefix'] + '_remixt.h5'
    remixt_raw_dir = args['output_prefix'] + '_raw_dir'

    remixt_brk_cn_csv = args['output_prefix'] + '_remixt_brk_cn.csv.gz'
    remixt_cn_csv = args['output_prefix'] + '_remixt_cn.csv.gz'
    remixt_minor_modes_csv = args['output_prefix'] + '_remixt_minor_modes.csv.gz'
    remixt_mix_csv = args['output_prefix'] + '_remixt_mix.csv.gz'
    remixt_read_depth_csv = args['output_prefix'] + '_remixt_read_depth.csv.gz'
    remixt_stats_csv = args['output_prefix'] + '_remixt_stats.csv.gz'

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    if run_remixt:
        breakpoints = inputs['breakpoints']
        workflow.subworkflow(
            name='remixt',
            func=remixt.create_remixt_workflow,
            args=(
                mgd.InputFile(tumour, extensions=['.bai']),
                mgd.InputFile(normal, extensions=['.bai']),
                mgd.InputFile(breakpoints),
                sample_id,
                mgd.OutputFile(remixt_outfile),
                mgd.OutputFile(remixt_brk_cn_csv),
                mgd.OutputFile(remixt_cn_csv),
                mgd.OutputFile(remixt_minor_modes_csv),
                mgd.OutputFile(remixt_mix_csv),
                mgd.OutputFile(remixt_read_depth_csv),
                mgd.OutputFile(remixt_stats_csv),
                refdir_paths['refdata_remixt'],
                remixt_raw_dir,
                refdir_paths['reference'],
                chromosomes
            ),
            kwargs={'single_node': args['single_node'], 'sex': sex}
        )

    if run_titan:
        workflow.subworkflow(
            name='titan',
            func=titan.create_titan_workflow,
            axes=('sample_id',),
            args=(
                mgd.InputFile(tumour,extensions=['.bai']),
                mgd.InputFile(normal,extensions=['.bai']),
                mgd.InputFile(targets),
                mgd.OutputFile(titan_outfile),
                mgd.OutputFile(titan_params),
                mgd.OutputFile(titan_segs),
                mgd.OutputFile(titan_igv_segs),
                mgd.OutputFile(titan_parsed),
                mgd.OutputFile(titan_plots),
                mgd.OutputFile(titan_tar_outputs),
                mgd.OutputFile(museq_vcf),
                sample_id,
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

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'copynumber_calling'}
        }
    )

    pyp.run(workflow)
