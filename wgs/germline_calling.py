import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import germline_calling


def germline_calling_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    normal = inputs['normal']
    normal_id = inputs['normal_id']

    museq_ss_vcf = args['output_prefix'] + '_museq_single_annotated.vcf.gz'
    museq_ss_maf = args['output_prefix'] + '_museq_single_annotated.maf'
    museq_single_pdf = args['output_prefix'] + '_single_museqportrait.pdf'

    samtools_germline_vcf = args['output_prefix'] + '_samtools_germline.vcf.gz'
    samtools_germline_maf = args['output_prefix'] + '_samtools_germline.maf'
    samtools_roh = args['output_prefix'] + '_roh.csv.gz'

    freebayes_germline_vcf = args['output_prefix'] + '_freebayes_germline.vcf.gz'
    freebayes_germline_maf = args['output_prefix'] + '_freebayes_germline.maf'

    rtg_germline_vcf = args['output_prefix'] + '_rtg_germline.vcf.gz'
    rtg_germline_maf = args['output_prefix'] + '_rtg_germline.maf'

    consensus_germline_maf = args['output_prefix'] + '_consensus_germline.maf'

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx()
    )

    workflow.subworkflow(
        name='germline_calling',
        func=germline_calling.create_germline_calling_workflow,
        args=(
            mgd.InputFile(normal, extensions=['.bai']),
            mgd.OutputFile(museq_ss_vcf),
            mgd.OutputFile(museq_ss_maf),
            mgd.OutputFile(museq_single_pdf),
            mgd.OutputFile(samtools_germline_vcf),
            mgd.OutputFile(samtools_germline_maf),
            mgd.OutputFile(samtools_roh),
            mgd.OutputFile(freebayes_germline_vcf),
            mgd.OutputFile(freebayes_germline_maf),
            mgd.OutputFile(rtg_germline_vcf),
            mgd.OutputFile(rtg_germline_maf),
            mgd.OutputFile(consensus_germline_maf),
            args['refdir'],
            normal_id
        ),
        kwargs={
            'single_node': args['single_node'],
        }
    )

    filenames = [
        museq_ss_vcf,
        museq_ss_maf,
        museq_single_pdf,
        samtools_germline_vcf,
        samtools_germline_maf,
        samtools_roh,
        freebayes_germline_vcf,
        freebayes_germline_maf,
        rtg_germline_vcf,
        rtg_germline_maf,
        consensus_germline_maf
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
