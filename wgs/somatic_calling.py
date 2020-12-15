import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import somatic_calling


def somatic_calling_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = list(tumours.keys())

    tumour_ids = helpers.get_values_from_input(inputs, 'tumour_id')
    normal_ids = helpers.get_values_from_input(inputs, 'normal_id')

    var_dir = os.path.join(args['out_dir'], 'somatic')
    museq_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_paired_annotated.vcf.gz')
    museq_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_paired_annotated.maf')
    museq_paired_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_paired_museqportrait.pdf')

    strelka_snv_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_snv_annotated.vcf.gz')
    strelka_snv_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_snv_annotated.maf')
    strelka_indel_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_indel_annotated.vcf.gz')
    strelka_indel_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_strelka_indel_annotated.maf')

    mutect_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_mutect.vcf.gz')
    mutect_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_mutect.maf')

    consensus_somatic_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_somatic.maf')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    workflow.subworkflow(
        name='variant_calling',
        func=somatic_calling.create_somatic_calling_workflow,
        args=(
            samples,
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('museq_vcf', 'sample_id', template=museq_vcf, axes_origin=[]),
            mgd.OutputFile('museq_maf', 'sample_id', template=museq_maf, axes_origin=[]),
            mgd.OutputFile('museq_paired_pdf', 'sample_id', template=museq_paired_pdf, axes_origin=[]),
            mgd.OutputFile('strelka_snv_vcf', 'sample_id', template=strelka_snv_vcf, axes_origin=[]),
            mgd.OutputFile('strelka_snv_maf', 'sample_id', template=strelka_snv_maf, axes_origin=[]),
            mgd.OutputFile('strelka_indel_vcf', 'sample_id', template=strelka_indel_vcf, axes_origin=[]),
            mgd.OutputFile('strelka_indel_maf', 'sample_id', template=strelka_indel_maf, axes_origin=[]),
            mgd.OutputFile('mutect_vcf', 'sample_id', template=mutect_vcf, axes_origin=[]),
            mgd.OutputFile('mutect_maf', 'sample_id', template=mutect_maf, axes_origin=[]),
            mgd.OutputFile('consensus_somatic_maf', 'sample_id', template=consensus_somatic_maf, axes_origin=[]),
            args['refdir'],
            normal_ids,
            tumour_ids,
        ),
        kwargs={
            'single_node': args['single_node'],
            'is_exome': args['is_exome'],
        }
    )

    filenames = [
        museq_vcf, museq_maf, museq_paired_pdf,
        strelka_snv_vcf, strelka_snv_maf,
        strelka_indel_vcf, strelka_indel_maf,
        mutect_vcf, mutect_maf,
        consensus_somatic_maf
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
