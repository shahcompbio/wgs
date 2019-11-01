import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import remixt
from wgs.workflows import titan

from alignment import paired_alignment
from sv_calling import call_breakpoints
from variant_calling import call_variants


def wgs_workflow(args):
    run_var_calling = args['variant_calling']
    run_cn_calling = args['copynumber_calling']
    run_bkp_calling = args['breakpoint_calling']
    run_alignment = args['alignment']

    if not any((run_var_calling, run_cn_calling, run_bkp_calling, run_alignment)):
        run_cn_calling = True
        run_bkp_calling = True
        run_var_calling = True
        run_alignment = True

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    targets = helpers.get_values_from_input(inputs, 'target_list')
    samples = tumours.keys()

    single_node = args['single_node']

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    if run_alignment:
        tumour_fastqs_r1, tumour_fastqs_r2 = helpers.get_fastqs(inputs, samples, 'tumour')
        normal_fastqs_r1, normal_fastqs_r2 = helpers.get_fastqs(inputs, samples, 'normal')

        workflow.subworkflow(
            name='wgs_alignment_paired_lanes',
            func=paired_alignment,
            ctx={'docker_image': config['cna_calling']['docker']['wgs']},
            args=(
                config,
                mgd.OutputFile("tumour.bam", 'sample_id', fnames=tumours,
                               extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile("normal.bam", 'sample_id', fnames=normals,
                               extensions=['.bai'], axes_origin=[]),
                samples,
                tumour_fastqs_r1,
                tumour_fastqs_r2,
                normal_fastqs_r1,
                normal_fastqs_r2,
                args['out_dir'],
            ),
            kwargs={'single_node': single_node}
        )

    if run_var_calling:
        museq_dir = os.path.join(args['out_dir'], 'variants')
        museq_vcf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_museq_paired_annotated.vcf.gz')
        museq_ss_vcf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_museq_single_annotated.vcf.gz')
        strelka_snv_vcf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_strelka_snv_annotated.vcf.gz')
        strelka_indel_vcf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_strelka_indel_annotated.vcf.gz')
        parsed_snv_csv = os.path.join(museq_dir, '{sample_id}', '{sample_id}_allcalls.csv')
        museq_paired_pdf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_paired_museqportrait.pdf')
        museq_single_pdf = os.path.join(museq_dir, '{sample_id}', '{sample_id}_single_museqportrait.pdf')
        workflow.subworkflow(
            name='variant_calling',
            func=call_variants,
            ctx={'docker_image': config['cna_calling']['docker']['wgs']},
            args=(
                samples,
                config,
                mgd.OutputFile('parsed_snv_csv', 'sample_id', template=parsed_snv_csv, axes_origin=[]),
                mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                              extensions=['.bai'], axes_origin=[]),
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile('museq', 'sample_id', template=museq_vcf, axes_origin=[]),
                mgd.OutputFile('museq_ss', 'sample_id', template=museq_ss_vcf, axes_origin=[]),
                mgd.OutputFile('strelka_snv', 'sample_id', template=strelka_snv_vcf, axes_origin=[]),
                mgd.OutputFile('strelka_indel', 'sample_id', template=strelka_indel_vcf, axes_origin=[]),
                mgd.OutputFile('museq_paired_pdf', 'sample_id', template=museq_paired_pdf, axes_origin=[]),
                mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
            ),
            kwargs={'single_node': single_node}
        )

    if run_bkp_calling:
        sv_outdir = os.path.join(args['out_dir'], 'breakpoints', '{sample_id}')
        destruct_breakpoints = os.path.join(sv_outdir, '{sample_id}_destruct_breakpoints.csv')
        destruct_library = os.path.join(sv_outdir, '{sample_id}_destruct_library.csv')
        destruct_raw_breakpoints = os.path.join(sv_outdir, '{sample_id}_destruct_raw_breakpoints.csv')
        destruct_raw_library = os.path.join(sv_outdir, '{sample_id}_destruct_raw_library.csv')
        destruct_reads = os.path.join(sv_outdir, '{sample_id}_destruct_reads.csv')
        lumpy_vcf = os.path.join(sv_outdir, '{sample_id}_lumpy.vcf')
        parsed_csv = os.path.join(sv_outdir, '{sample_id}_filtered_consensus_calls.csv')
        workflow.subworkflow(
            name="call_breakpoints",
            func=call_breakpoints,
            ctx={'docker_image': config['sv_calling']['docker']['wgs']},
            args=(
                samples,
                config,
                mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                              extensions=['.bai'], axes_origin=[]),
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile('destruct_raw_breakpoints', 'sample_id', template=destruct_raw_breakpoints,
                               axes_origin=[]),
                mgd.OutputFile('destruct_raw_library', 'sample_id', template=destruct_raw_library, axes_origin=[]),
                mgd.OutputFile('destruct_breakpoints', 'sample_id', template=destruct_breakpoints, axes_origin=[]),
                mgd.OutputFile('destruct_library', 'sample_id', template=destruct_library, axes_origin=[]),
                mgd.OutputFile('destruct_reads', 'sample_id', template=destruct_reads, axes_origin=[]),
                mgd.OutputFile('lumpy_vcf', 'sample_id', template=lumpy_vcf, axes_origin=[]),
                mgd.OutputFile('parsed_csv', 'sample_id', template=parsed_csv, axes_origin=[])
            ),
            kwargs={'single_node': single_node}
        )

    if run_cn_calling:
        if run_bkp_calling:
            sv_outdir = os.path.join(args['out_dir'], 'breakpoints', '{sample_id}')
            destruct_breakpoints = os.path.join(sv_outdir, 'destruct_breakpoints.csv')
            destruct_breakpoints = mgd.InputFile(
                'destruct_breakpoints', 'sample_id', axes_origin=[],
                template=destruct_breakpoints
            )
        else:
            destruct_breakpoints = None

        cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')
        remixt_raw_dir = os.path.join(cna_outdir, 'remixt', 'raw_data')
        titan_raw_dir = os.path.join(cna_outdir, 'titan')
        remixt_results_filename = os.path.join(cna_outdir, 'remixt', '{sample_id}_results.h5')
        titan_segments_filename = os.path.join(titan_raw_dir, '{sample_id}_segments.h5')
        titan_markers_filename = os.path.join(titan_raw_dir, '{sample_id}_markers.h5')
        titan_params_filename = os.path.join(titan_raw_dir, '{sample_id}_params.h5')
        workflow.subworkflow(
            name='titan',
            func=titan.create_titan_workflow,
            ctx={'docker_image': config['cna_calling']['docker']['wgs']},
            axes=('sample_id',),
            args=(
                mgd.InputFile('tumour.bam', 'sample_id', fnames=tumours, extensions=['.bai']),
                mgd.InputFile('normal.bam', 'sample_id', fnames=normals, extensions=['.bai']),
                mgd.InputFile("target_list", 'sample_id', fnames=targets, axes_origin=[]),
                mgd.Template(titan_raw_dir, 'sample_id'),
                mgd.OutputFile('titan_segments_filename', 'sample_id', axes_origin=[],
                               template=titan_segments_filename),
                mgd.OutputFile('titan_params_filename', 'sample_id', axes_origin=[], template=titan_params_filename),
                mgd.OutputFile('titan_markers_filename', 'sample_id', axes_origin=[], template=titan_markers_filename),
                config['globals'],
                config['cna_calling'],
                config['cna_calling']['titan_intervals'],
                mgd.InputInstance('sample_id'),
            ),
            kwargs={'single_node': single_node}
        )
        workflow.subworkflow(
            name='remixt',
            func=remixt.create_remixt_workflow,
            ctx={'docker_image': config['cna_calling']['docker']['wgs']},
            axes=('sample_id',),
            args=(
                mgd.InputFile('tumour.bam', 'sample_id', fnames=tumours, extensions=['.bai']),
                mgd.InputFile('normal.bam', 'sample_id', fnames=normals, extensions=['.bai']),
                destruct_breakpoints,
                mgd.InputInstance('sample_id'),
                config['cna_calling']['remixt_refdata'],
                mgd.OutputFile('remixt_results_filename', 'sample_id', axes_origin=[],
                               template=remixt_results_filename),
                mgd.Template(remixt_raw_dir, 'sample_id'),
                config['cna_calling']['min_num_reads'],
                config['globals']
            ),
            kwargs={'docker_containers': config['cna_calling']['docker']}
        )

    pyp.run(workflow)
