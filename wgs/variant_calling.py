import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def call_germlines_only(
        samples, normals, museq_ss_vcf, samtools_germline_vcf, roh_calls,
        museq_single_pdf, refdir, single_node=False
):
    museq_ss_vcf = dict([(sampid, museq_ss_vcf[sampid])
                         for sampid in samples])
    museq_single_pdf = dict([(sampid, museq_single_pdf[sampid])
                             for sampid in samples])
    samtools_germline_vcf = dict([(sampid, samtools_germline_vcf[sampid])
                                  for sampid in samples])
    roh_calls = dict([(sampid, roh_calls[sampid])
                      for sampid in samples])

    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    paths_refdir = config.refdir_data(refdir)['paths']

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
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
            paths_refdir['reference'],
            chromosomes,
        ),
        kwargs={
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
            'germline_refdata': paths_refdir['germline_portrait_ref'],
            'thousand_genomes': paths_refdir['thousand_genomes'],
            'dbsnp': paths_refdir['dbsnp'],
        }
    )

    workflow.subworkflow(
        name="samtools_germline",
        func='wgs.workflows.samtools_germline.create_samtools_germline_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("samtools_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile("roh_calls.csv", 'sample_id',
                           fnames=roh_calls),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            paths_refdir['reference'],
            chromosomes
        ),
        kwargs={
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
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
                }
    )

    workflow.subworkflow(
        name="annotate_germline_samtools",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("samtools_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile("samtools_germlines_anno.vcf.gz", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=samtools_germline_vcf),
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
                }
    )


    return workflow


def call_variants(
        samples,
        somatic_calls, somatic_snpeff, somatic_ma, somatic_ids,
        indel_calls, indel_snpeff, indel_ma, indel_ids,
        germline_calls, germline_snpeff, germline_ma, germline_ids,
        tumours, normals, museq_vcf, museq_ss_vcf, samtools_germlines_vcf, roh_calls,
        strelka_snv_vcf, strelka_indel_vcf,
        museq_paired_pdf, museq_single_pdf, refdir,
        single_node=False, strelka_depth_filter=True
):
    strelka_snv_vcf = dict([(sampid, strelka_snv_vcf[sampid])
                            for sampid in samples])
    strelka_indel_vcf = dict([(sampid, strelka_indel_vcf[sampid])
                              for sampid in samples])
    museq_vcf = dict([(sampid, museq_vcf[sampid])
                      for sampid in samples])
    museq_ss_vcf = dict([(sampid, museq_ss_vcf[sampid])
                         for sampid in samples])
    samtools_germlines_vcf = dict([(sampid, samtools_germlines_vcf[sampid])
                                   for sampid in samples])
    roh_calls = dict([(sampid, roh_calls[sampid])
                      for sampid in samples])

    museq_paired_pdf = dict([(sampid, museq_paired_pdf[sampid])
                             for sampid in samples])
    museq_single_pdf = dict([(sampid, museq_single_pdf[sampid])
                             for sampid in samples])

    somatic_calls = dict([(sampid, somatic_calls[sampid])
                          for sampid in samples])
    somatic_snpeff = dict([(sampid, somatic_snpeff[sampid])
                           for sampid in samples])
    somatic_ma = dict([(sampid, somatic_ma[sampid])
                       for sampid in samples])
    somatic_ids = dict([(sampid, somatic_ids[sampid])
                        for sampid in samples])

    indel_calls = dict([(sampid, indel_calls[sampid])
                        for sampid in samples])
    indel_snpeff = dict([(sampid, indel_snpeff[sampid])
                         for sampid in samples])
    indel_ma = dict([(sampid, indel_ma[sampid])
                     for sampid in samples])
    indel_ids = dict([(sampid, indel_ids[sampid])
                      for sampid in samples])

    germline_calls = dict([(sampid, germline_calls[sampid])
                           for sampid in samples])
    germline_snpeff = dict([(sampid, germline_snpeff[sampid])
                            for sampid in samples])
    germline_ma = dict([(sampid, germline_ma[sampid])
                        for sampid in samples])
    germline_ids = dict([(sampid, germline_ids[sampid])
                         for sampid in samples])

    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    paths_refdir = config.refdir_data(refdir)['paths']

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
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
            paths_refdir['reference'],
            chromosomes
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
            paths_refdir['reference'],
            chromosomes
        ),
        kwargs={
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
            'germline_refdata': paths_refdir['germline_portrait_ref'],
            'thousand_genomes': paths_refdir['thousand_genomes'],
            'dbsnp': paths_refdir['dbsnp'],
        }
    )

    workflow.subworkflow(
        name="samtools_germline",
        func='wgs.workflows.samtools_germline.create_samtools_germline_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("samtools_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile("roh_calls.csv.gz", 'sample_id',
                           fnames=roh_calls),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            paths_refdir['reference'],
            chromosomes
        ),
        kwargs={
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
            mgd.TempOutputFile('strelka_indel.vcf.gz', 'sample_id'),
            mgd.TempOutputFile('strelka_snv.vcf.gz', 'sample_id'),
            paths_refdir['reference'],
            chromosomes
        ),
        kwargs={
            'single_node': single_node,
            'use_depth_thresholds': strelka_depth_filter
        },

    )

    workflow.subworkflow(
        name="annotate_paired_museq",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_vcf),
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
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
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
                }
    )


    workflow.subworkflow(
        name="annotate_germline_samtools",
        func='wgs.workflows.vcf_annotation.create_annotation_workflow',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("samtools_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile("samtools_germlines_ann.vcf.gz", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=samtools_germlines_vcf),
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
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
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
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
            paths_refdir['snpeff_config'],
            paths_refdir['mutation_assessor'],
            paths_refdir['dbsnp'],
            paths_refdir['thousand_genomes'],
            paths_refdir['cosmic'],
            paths_refdir['blacklist']
        ),
        kwargs={'vcftools_docker': config.containers('vcftools'),
                'snpeff_docker': config.containers('vcftools'),
                }
    )

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
            mgd.OutputFile('somatic_snpeff', 'sample_id', fnames=somatic_snpeff),
            mgd.OutputFile('somatic_ma', 'sample_id', fnames=somatic_ma),
            mgd.OutputFile('somatic_ids', 'sample_id', fnames=somatic_ids),
            mgd.OutputFile('indel_csv', 'sample_id', fnames=indel_calls),
            mgd.OutputFile('indel_snpeff', 'sample_id', fnames=indel_snpeff),
            mgd.OutputFile('indel_ma', 'sample_id', fnames=indel_ma),
            mgd.OutputFile('indel_ids', 'sample_id', fnames=indel_ids),
            mgd.OutputFile('germline_csv', 'sample_id', fnames=germline_calls),
            mgd.OutputFile('germline_snpeff', 'sample_id', fnames=germline_snpeff),
            mgd.OutputFile('germline_ma', 'sample_id', fnames=germline_ma),
            mgd.OutputFile('germline_ids', 'sample_id', fnames=germline_ids),
            refdir,
        ),
    )

    return workflow


def variant_calling_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = tumours.keys()

    var_dir = os.path.join(args['out_dir'], 'variants')
    museq_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_paired_annotated.vcf.gz')
    museq_ss_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_single_annotated.vcf.gz')

    samtools_germline_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_samtools_germline.vcf.gz')
    samtools_roh = os.path.join(var_dir, '{sample_id}', '{sample_id}_roh.csv')

    strelka_snv_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_snv_annotated.vcf.gz')
    strelka_indel_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_indel_annotated.vcf.gz')
    museq_paired_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_paired_museqportrait.pdf')
    museq_single_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_single_museqportrait.pdf')

    somatic_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_somatic.csv.gz')
    somatic_snpeff = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_somatic_snpeff.csv.gz')
    somatic_ma = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_somatic_ma.csv.gz')
    somatic_ids = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_somatic_ids.csv.gz')

    indel_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_indel.csv.gz')
    indel_snpeff = os.path.join(var_dir, '{sample_id}', '{sample_id}_indel_snpeff.csv.gz')
    indel_ma = os.path.join(var_dir, '{sample_id}', '{sample_id}_indel_ma.csv.gz')
    indel_ids = os.path.join(var_dir, '{sample_id}', '{sample_id}_indel_ids.csv.gz')

    germline_csv = os.path.join(var_dir, '{sample_id}', '{sample_id}_germline.csv.gz')
    germline_snpeff = os.path.join(var_dir, '{sample_id}', '{sample_id}_germline_snpeff.csv.gz')
    germline_ma = os.path.join(var_dir, '{sample_id}', '{sample_id}_germline_ma.csv.gz')
    germline_ids = os.path.join(var_dir, '{sample_id}', '{sample_id}_germline_ids.csv.gz')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
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
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile('museq_ss', 'sample_id', template=museq_ss_vcf, axes_origin=[]),
                mgd.OutputFile('samtools_germline', 'sample_id', template=samtools_germline_vcf, axes_origin=[]),
                mgd.OutputFile('samtools_roh', 'sample_id', template=samtools_roh, axes_origin=[]),
                mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
                args['refdir']
            ),
            kwargs={'single_node': args['single_node']}
        )
    else:
        workflow.subworkflow(
            name='variant_calling',
            func=call_variants,
            args=(
                samples,
                mgd.OutputFile('somatic_csv', 'sample_id', template=somatic_csv, axes_origin=[]),
                mgd.OutputFile('somatic_snpeff', 'sample_id', template=somatic_snpeff, axes_origin=[]),
                mgd.OutputFile('somatic_ma', 'sample_id', template=somatic_ma, axes_origin=[]),
                mgd.OutputFile('somatic_ids', 'sample_id', template=somatic_ids, axes_origin=[]),
                mgd.OutputFile('indel_csv', 'sample_id', template=indel_csv, axes_origin=[]),
                mgd.OutputFile('indel_snpeff', 'sample_id', template=indel_snpeff, axes_origin=[]),
                mgd.OutputFile('indel_ma', 'sample_id', template=indel_ma, axes_origin=[]),
                mgd.OutputFile('indel_ids', 'sample_id', template=indel_ids, axes_origin=[]),
                mgd.OutputFile('germline_csv', 'sample_id', template=germline_csv, axes_origin=[]),
                mgd.OutputFile('germline_snpeff', 'sample_id', template=germline_snpeff, axes_origin=[]),
                mgd.OutputFile('germline_ma', 'sample_id', template=germline_ma, axes_origin=[]),
                mgd.OutputFile('germline_ids', 'sample_id', template=germline_ids, axes_origin=[]),
                mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                              extensions=['.bai'], axes_origin=[]),
                mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                              extensions=['.bai'], axes_origin=[]),
                mgd.OutputFile('museq', 'sample_id', template=museq_vcf, axes_origin=[]),
                mgd.OutputFile('museq_ss', 'sample_id', template=museq_ss_vcf, axes_origin=[]),
                mgd.OutputFile('samtools_germline', 'sample_id', template=samtools_germline_vcf, axes_origin=[]),
                mgd.OutputFile('roh_calls', 'sample_id', template=samtools_roh, axes_origin=[]),
                mgd.OutputFile('strelka_snv', 'sample_id', template=strelka_snv_vcf, axes_origin=[]),
                mgd.OutputFile('strelka_indel', 'sample_id', template=strelka_indel_vcf, axes_origin=[]),
                mgd.OutputFile('museq_paired_pdf', 'sample_id', template=museq_paired_pdf, axes_origin=[]),
                mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
                args['refdir'],
            ),
            kwargs={
                'single_node': args['single_node'],
                'strelka_depth_filter': not args['remove_strelka_depth_filter'],
            }
        )

        filenames = [
            somatic_csv, somatic_snpeff, somatic_ma, somatic_ids,
            indel_csv, indel_snpeff, indel_ma, indel_ids,
            germline_csv, germline_snpeff, germline_ma, germline_ids,
            museq_vcf, museq_ss_vcf, strelka_snv_vcf, strelka_indel_vcf,
            museq_paired_pdf, museq_single_pdf
        ]

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
