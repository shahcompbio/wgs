from variant_calling import call_variants
from sv_calling import call_breakpoints
from cna_calling import call_copynumber
import pypeliner
from wgs.utils import helpers
import os
import pypeliner.managed as mgd
from alignment import paired_alignment
from workflows import titan
from workflows import remixt

def get_fastqs(inputs, samples, sample_type):
    fq1 = {}
    fq2 = {}

    for sample in samples:
        fastqs = inputs[sample]['fastqs'][sample_type]
        for lane in fastqs:
            fq1[(sample, lane)] = fastqs[lane]['fastq1']
            fq2[(sample, lane)] = fastqs[lane]['fastq2']

    return fq1, fq2


def wgs_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    if args['alignment']:
        tumour_fastqs_r1, tumour_fastqs_r2 = get_fastqs(inputs, samples, 'tumour')
        normal_fastqs_r1, normal_fastqs_r2 = get_fastqs(inputs, samples, 'normal')

        normal_alignment_template = os.path.join(
            args['out_dir'], 'alignment', '{norm_sample_id}', '{norm_lane}', 'normal'
        )
        tumour_alignment_template = os.path.join(
            args['out_dir'], 'alignment', '{tum_sample_id}', '{tum_lane}', 'tumour'
        )

        workflow.subworkflow(
            name='wgs_alignment_paired_lanes',
            func=paired_alignment,
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
                normal_alignment_template,
                tumour_alignment_template,
            )
        )

    museq_dir = os.path.join(args['out_dir'], 'variants')
    museq_vcf = os.path.join(museq_dir, '{sample_id}', 'museq_paired_annotated.vcf.gz')
    museq_ss_vcf = os.path.join(museq_dir, '{sample_id}', 'museq_single_annotated.vcf.gz')
    strelka_snv_vcf = os.path.join(museq_dir, '{sample_id}', 'strelka_snv_annotated.vcf.gz')
    strelka_indel_vcf = os.path.join(museq_dir, '{sample_id}', 'strelka_indel_annotated.vcf.gz')
    parsed_snv_csv = os.path.join(museq_dir, '{sample_id}', 'allcalls.csv')
    museq_paired_pdf = os.path.join(museq_dir, '{sample_id}', 'paired_museqportrait.pdf')
    museq_paired_pdf_txt = os.path.join(museq_dir, '{sample_id}', 'paired_museqportrait.txt')
    museq_single_pdf = os.path.join(museq_dir, '{sample_id}', 'single_museqportrait.pdf')
    museq_single_pdf_txt = os.path.join(museq_dir, '{sample_id}', 'single_museqportrait.txt')
    workflow.subworkflow(
        name='variant_calling',
        func=call_variants,
        args=(
            samples,
            museq_dir,
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
            mgd.OutputFile('museq_paired_pdf_txt', 'sample_id', template=museq_paired_pdf_txt, axes_origin=[]),
            mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
            mgd.OutputFile('museq_single_pdf_txt', 'sample_id', template=museq_single_pdf_txt, axes_origin=[]),
        )
    )

    sv_outdir = os.path.join(args['out_dir'], 'breakpoints', '{sample_id}')
    destruct_breakpoints = os.path.join(sv_outdir, 'destruct_breakpoints.csv')
    destruct_library = os.path.join(sv_outdir, 'destruct_library.csv')
    destruct_raw_breakpoints = os.path.join(sv_outdir, 'destruct_raw_breakpoints.csv')
    destruct_raw_library = os.path.join(sv_outdir, 'destruct_raw_library.csv')
    destruct_reads = os.path.join(sv_outdir, 'destruct_reads.csv')
    lumpy_vcf = os.path.join(sv_outdir, 'lumpy.vcf')
    parsed_csv = os.path.join(sv_outdir, 'filtered_consensus_calls.csv')
    workflow.subworkflow(
        name="call_breakpoints",
        func=call_breakpoints,
        args=(
            samples,
            config,
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('destruct_raw_breakpoints', 'sample_id', template=destruct_raw_breakpoints, axes_origin=[]),
            mgd.OutputFile('destruct_raw_library', 'sample_id', template=destruct_raw_library, axes_origin=[]),
            mgd.OutputFile('destruct_breakpoints', 'sample_id', template=destruct_breakpoints, axes_origin=[]),
            mgd.OutputFile('destruct_library', 'sample_id', template=destruct_library, axes_origin=[]),
            mgd.OutputFile('destruct_reads', 'sample_id', template=destruct_reads, axes_origin=[]),
            mgd.OutputFile('lumpy_vcf', 'sample_id', template=lumpy_vcf, axes_origin=[]),
            mgd.OutputFile('parsed_csv', 'sample_id', template=parsed_csv, axes_origin=[])
        )
    )

    cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')
    remixt_raw_dir = os.path.join(cna_outdir, 'remixt', 'raw_data')
    titan_raw_dir = os.path.join(cna_outdir, 'titan')
    remixt_results_filename = os.path.join(cna_outdir, 'remixt', 'results.h5')
    titan_segments_filename = os.path.join(titan_raw_dir, 'segments.h5')
    titan_markers_filename = os.path.join(titan_raw_dir, 'markers.h5')
    titan_params_filename = os.path.join(titan_raw_dir, 'params.h5')
    workflow.subworkflow(
        name='titan',
        func=titan.create_titan_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.Template(titan_raw_dir, 'sample_id'),
            mgd.OutputFile('titan_segments_filename', 'sample_id', axes_origin=[], template=titan_segments_filename),
            mgd.OutputFile('titan_params_filename', 'sample_id', axes_origin=[], template=titan_params_filename),
            mgd.OutputFile('titan_markers_filename', 'sample_id', axes_origin=[], template=titan_markers_filename),
            config['globals'],
            config['cna_calling'],
            config['cna_calling']['titan_intervals'],
            mgd.InputInstance('sample_id'),
        ),
    )
    workflow.subworkflow(
        name='remixt',
        func=remixt.create_remixt_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('destruct_breakpoints', 'sample_id', axes_origin=[], template=destruct_breakpoints),
            mgd.InputInstance('sample_id'),
            config['cna_calling']['remixt_refdata'],
            mgd.OutputFile('remixt_results_filename', 'sample_id', axes_origin=[], template=remixt_results_filename),
            mgd.Template(remixt_raw_dir, 'sample_id'),
            config['cna_calling']['min_num_reads']
        ),
    )

    pyp.run(workflow)