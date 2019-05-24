import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import mutationseq
from wgs.workflows import strelka
from wgs.workflows import variant_calling_consensus
from wgs.workflows import vcf_annotation


def call_variants(
        samples, config, parsed_csv,
        tumours, normals, museq_vcf, museq_ss_vcf,
        strelka_snv_vcf, strelka_indel_vcf,
        museq_paired_pdf, museq_single_pdf,
        single_node=False
):
    strelka_snv_vcf = dict([(sampid, strelka_snv_vcf[sampid])
                            for sampid in samples])
    strelka_indel_vcf = dict([(sampid, strelka_indel_vcf[sampid])
                              for sampid in samples])
    museq_vcf = dict([(sampid, museq_vcf[sampid])
                      for sampid in samples])
    museq_ss_vcf = dict([(sampid, museq_ss_vcf[sampid])
                         for sampid in samples])
    museq_paired_pdf = dict([(sampid, museq_paired_pdf[sampid])
                             for sampid in samples])
    museq_single_pdf = dict([(sampid, museq_single_pdf[sampid])
                             for sampid in samples])
    parsed_csv = dict([(sampid, parsed_csv[sampid])
                       for sampid in samples])

    ctx = {'docker_image': config['variant_calling']['docker']['wgs']}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)


    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="mutationseq_paired",
        func=mutationseq.create_museq_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_paired_pdf', 'sample_id', fnames=museq_paired_pdf),
            config['globals'],
            config['variant_calling'],
        ),
        kwargs={
            'tumour_bam': mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                                        extensions=['.bai'], axes_origin=[]),
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="mutationseq_single",
        func=mutationseq.create_museq_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_single_pdf', 'sample_id', fnames=museq_single_pdf),
            config['globals'],
            config['variant_calling'],
        ),
        kwargs={
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="strelka",
        func=strelka.create_strelka_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            config['reference'],
            mgd.TempOutputFile('strelka_indel.vcf.gz', 'sample_id'),
            mgd.TempOutputFile('strelka_snv.vcf.gz', 'sample_id'),
            config['globals'],
            config['variant_calling'],
        ),
        kwargs={'single_node': single_node},
    )

    workflow.subworkflow(
        name="annotate_paired_museq",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['variant_calling']['docker']['vcftools'],
                'snpeff_docker': config['variant_calling']['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_germline_museq",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_germlines_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_ss_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['variant_calling']['docker']['vcftools'],
                'snpeff_docker': config['variant_calling']['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_strelka",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_snv_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['variant_calling']['docker']['vcftools'],
                'snpeff_docker': config['variant_calling']['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_strelka_indel",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_indel.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_indel_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_indel_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['variant_calling']['docker']['vcftools'],
                'snpeff_docker': config['variant_calling']['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="consensus_calling",
        func=variant_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("museq_germlines_ann.vcf.gz", 'sample_id', fnames=museq_ss_vcf),
            mgd.InputFile("museq_snv_ann.vcf.gz", 'sample_id', fnames=museq_vcf),
            mgd.InputFile("strelka_snv_ann.vcf.gz", 'sample_id', fnames=strelka_snv_vcf),
            mgd.InputFile("strelka_indel_ann.vcf.gz", 'sample_id', fnames=strelka_indel_vcf),
            mgd.OutputFile('parsed_csv', 'sample_id', fnames=parsed_csv),
            config['globals'],
            config['variant_calling'],
        ),
    )

    return workflow


def variant_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = tumours.keys()

    museq_dir = os.path.join(args['out_dir'], 'variants')
    museq_vcf = os.path.join(museq_dir, '{sample_id}', 'museq_paired_annotated.vcf.gz')
    museq_ss_vcf = os.path.join(museq_dir, '{sample_id}', 'museq_single_annotated.vcf.gz')
    strelka_snv_vcf = os.path.join(museq_dir, '{sample_id}', 'strelka_snv_annotated.vcf.gz')
    strelka_indel_vcf = os.path.join(museq_dir, '{sample_id}', 'strelka_indel_annotated.vcf.gz')
    parsed_snv_csv = os.path.join(museq_dir, '{sample_id}', 'allcalls.csv')
    museq_paired_pdf = os.path.join(museq_dir, '{sample_id}', 'paired_museqportrait.pdf')
    museq_single_pdf = os.path.join(museq_dir, '{sample_id}', 'single_museqportrait.pdf')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    workflow.subworkflow(
        name='variant_calling',
        func=call_variants,
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
        kwargs={'single_node': args['single_node']}
    )

    pyp.run(workflow)
