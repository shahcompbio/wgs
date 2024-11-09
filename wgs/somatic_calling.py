import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import somatic_calling


def somatic_calling_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    tumour = inputs['tumour']
    normal = inputs['normal']
    tumour_id = inputs['tumour_id']
    normal_id = inputs['normal_id']

    museq_vcf = args['output_prefix'] + '_museq_paired_annotated.vcf.gz'
    museq_maf = args['output_prefix'] + '_museq_paired_annotated.maf'
    museq_paired_pdf = args['output_prefix'] + '_paired_museqportrait.pdf'

    strelka_snv_vcf = args['output_prefix'] + '_strelka_snv_annotated.vcf.gz'
    strelka_snv_maf = args['output_prefix'] + '_strelka_snv_annotated.maf'
    strelka_indel_vcf = args['output_prefix'] + '_strelka_indel_annotated.vcf.gz'
    strelka_indel_maf = args['output_prefix'] + '_strelka_indel_annotated.maf'

    mutect_vcf = args['output_prefix'] + '_mutect.vcf.gz'
    mutect_maf = args['output_prefix'] + '_mutect.maf'

    consensus_somatic_maf = args['output_prefix'] + '_consensus_somatic.maf'

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='variant_calling',
        func=somatic_calling.create_somatic_calling_workflow,
        args=(
            mgd.InputFile(tumour, extensions=['.bai']),
            mgd.InputFile(normal, extensions=['.bai']),
            mgd.OutputFile(museq_vcf),
            mgd.OutputFile(museq_maf),
            mgd.OutputFile(museq_paired_pdf),
            mgd.OutputFile(strelka_snv_vcf),
            mgd.OutputFile(strelka_snv_maf),
            mgd.OutputFile(strelka_indel_vcf),
            mgd.OutputFile(strelka_indel_maf),
            mgd.OutputFile(mutect_vcf),
            mgd.OutputFile(mutect_maf),
            mgd.OutputFile(consensus_somatic_maf),
            args['refdir'],
            normal_id,
            tumour_id,
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

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'variant_calling'}
        }
    )

    pyp.run(workflow)
