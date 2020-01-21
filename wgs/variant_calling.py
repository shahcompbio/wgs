import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def call_germlines_only(
        samples, config, normals, museq_ss_vcf,
        museq_single_pdf, single_node=False
):
    global_config = config['globals']
    config = config['variant_calling']

    museq_ss_vcf = dict([(sampid, museq_ss_vcf[sampid])
                         for sampid in samples])
    museq_single_pdf = dict([(sampid, museq_single_pdf[sampid])
                             for sampid in samples])

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs'])
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="mutationseq_single",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_single_pdf', 'sample_id', fnames=museq_single_pdf),
            global_config,
            config,
        ),
        kwargs={
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="annotate_germline_museq",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_germlines_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_ss_vcf),
            global_config,
            config['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['docker']['vcftools'],
                'snpeff_docker': config['docker']['vcftools'],
                }
    )

    return workflow


def call_variants(
        samples, config, somatic_calls, indel_calls, germline_calls, outdir,
        tumours, normals, museq_vcf, museq_ss_vcf,
        strelka_snv_vcf, strelka_indel_vcf,
        museq_paired_pdf, museq_single_pdf,
        single_node=False
):
    global_config = config['globals']
    config = config['variant_calling']

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

    somatic_calls = dict([(sampid, somatic_calls[sampid])
                          for sampid in samples])
    indel_calls = dict([(sampid, indel_calls[sampid])
                        for sampid in samples])
    germline_calls = dict([(sampid, germline_calls[sampid])
                           for sampid in samples])

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs'])
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="mutationseq_paired",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_paired_pdf', 'sample_id', fnames=museq_paired_pdf),
            global_config,
            config,
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
        func='wgs.workflows.mutationseq.create_museq_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_single_pdf', 'sample_id', fnames=museq_single_pdf),
            global_config,
            config,
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
        func='wgs.workflows.strelka.create_strelka_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            config['reference'],
            mgd.TempOutputFile('strelka_indel.vcf.gz', 'sample_id'),
            mgd.TempOutputFile('strelka_snv.vcf.gz', 'sample_id'),
            global_config,
            config,
        ),
        kwargs={'single_node': single_node},
    )

    workflow.subworkflow(
        name="annotate_paired_museq",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_vcf),
            global_config,
            config['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['docker']['vcftools'],
                'snpeff_docker': config['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_germline_museq",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_germlines_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_ss_vcf),
            global_config,
            config['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['docker']['vcftools'],
                'snpeff_docker': config['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_strelka",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_snv_vcf),
            global_config,
            config['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['docker']['vcftools'],
                'snpeff_docker': config['docker']['vcftools'],
                }
    )

    workflow.subworkflow(
        name="annotate_strelka_indel",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_indel.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_indel_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_indel_vcf),
            global_config,
            config['annotation_params'],
        ),
        kwargs={'vcftools_docker': config['docker']['vcftools'],
                'snpeff_docker': config['docker']['vcftools'],
                }
    )

    outdir = os.path.join(outdir, '{sample_id}')
    workflow.subworkflow(
        name="consensus_calling",
        func='wgs.workflows.variant_calling_consensus.create_consensus_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile("museq_germlines_ann.vcf.gz", 'sample_id', fnames=museq_ss_vcf),
            mgd.InputFile("museq_snv_ann.vcf.gz", 'sample_id', fnames=museq_vcf),
            mgd.InputFile("strelka_snv_ann.vcf.gz", 'sample_id', fnames=strelka_snv_vcf),
            mgd.InputFile("strelka_indel_ann.vcf.gz", 'sample_id', fnames=strelka_indel_vcf),
            mgd.OutputFile('somatic_csv', 'sample_id', fnames=somatic_calls),
            mgd.OutputFile('indel_csv', 'sample_id', fnames=indel_calls),
            mgd.OutputFile('germline_csv', 'sample_id', fnames=germline_calls),
            mgd.Template('template_outdir', 'sample_id', template=outdir),
            global_config,
            config,
        ),
    )

    return workflow


def variant_calling_workflow(args):
    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])
    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = tumours.keys()

    var_dir = os.path.join(args['out_dir'], 'variants')
    museq_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_paired_annotated.vcf.gz')
    museq_ss_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_single_annotated.vcf.gz')
    strelka_snv_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_snv_annotated.vcf.gz')
    strelka_indel_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_indel_annotated.vcf.gz')
    museq_paired_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_paired_museqportrait.pdf')
    museq_single_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_single_museqportrait.pdf')

    somatic_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_somatic.csv.gz')
    indel_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_indel.csv.gz')
    germline_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_germline.csv.gz')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config['variant_calling']['docker']['wgs'])
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    if not all(tumours.values()):
        workflow.subworkflow(
            name='variant_calling',
            func=call_germlines_only,
            args=(
                samples,
                config,
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile('museq_ss', 'sample_id', template=museq_ss_vcf, axes_origin=[]),
                mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
            ),
            kwargs={'single_node': args['single_node']}
        )
    else:
        workflow.subworkflow(
            name='variant_calling',
            func=call_variants,
            args=(
                samples,
                config,
                mgd.OutputFile('somatic_csv', 'sample_id', template=somatic_csv, axes_origin=[]),
                mgd.OutputFile('indel_csv', 'sample_id', template=indel_csv, axes_origin=[]),
                mgd.OutputFile('germline_csv', 'sample_id', template=germline_csv, axes_origin=[]),
                var_dir,
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

        filenames = [somatic_csv, indel_csv, germline_csv, museq_vcf,
            museq_ss_vcf, strelka_snv_vcf, strelka_indel_vcf,
            museq_paired_pdf, museq_single_pdf]

        outputted_filenames = helpers.expand_list(filenames, samples, "sample_id")

        workflow.transform(
            name='generate_meta_files_results',
            func='wgs.utils.helpers.generate_and_upload_metadata',
            args=(
                sys.argv[0:],
                args['out_dir'],
                outputted_filenames,
                mgd.OutputFile(meta_yaml)
            ),
            kwargs={
                'input_yaml_data': helpers.load_yaml(args['input_yaml']),
                'input_yaml': mgd.OutputFile(input_yaml_blob),
                'metadata': {'type': 'variant_calling'}
            }
        )

    pyp.run(workflow)
