import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import germline_calling


def germline_calling_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = list(normals.keys())

    normal_ids = helpers.get_values_from_input(inputs, 'normal_id')

    var_dir = os.path.join(args['out_dir'], 'germline')

    museq_ss_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_single_annotated.vcf.gz')
    museq_ss_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_museq_single_annotated.maf')
    museq_single_pdf = os.path.join(var_dir, '{sample_id}', '{sample_id}_single_museqportrait.pdf')

    samtools_germline_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_samtools_germline.vcf.gz')
    samtools_germline_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_samtools_germline.maf')
    samtools_roh = os.path.join(var_dir, '{sample_id}', '{sample_id}_roh.csv.gz')

    freebayes_germline_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_freebayes_germline.vcf.gz')
    freebayes_germline_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_freebayes_germline.maf')

    rtg_germline_vcf = os.path.join(var_dir, '{sample_id}', '{sample_id}_rtg_germline.vcf.gz')
    rtg_germline_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_rtg_germline.maf')

    consensus_germline_maf = os.path.join(var_dir, '{sample_id}', '{sample_id}_consensus_germline.maf')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx()
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    workflow.subworkflow(
        name='germline_calling',
        func=germline_calling.create_germline_calling_workflow,
        args=(
            samples,
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('museq_ss_vcf', 'sample_id', template=museq_ss_vcf, axes_origin=[]),
            mgd.OutputFile('museq_ss_maf', 'sample_id', template=museq_ss_maf, axes_origin=[]),
            mgd.OutputFile('museq_single_pdf', 'sample_id', template=museq_single_pdf, axes_origin=[]),
            mgd.OutputFile('samtools_germline_vcf', 'sample_id', template=samtools_germline_vcf, axes_origin=[]),
            mgd.OutputFile('samtools_germline_maf', 'sample_id', template=samtools_germline_maf, axes_origin=[]),
            mgd.OutputFile('samtools_roh', 'sample_id', template=samtools_roh, axes_origin=[]),
            mgd.OutputFile('freebayes_germline_vcf', 'sample_id', template=freebayes_germline_vcf, axes_origin=[]),
            mgd.OutputFile('freebayes_germline_maf', 'sample_id', template=freebayes_germline_maf, axes_origin=[]),
            mgd.OutputFile('rtg_germline_vcf', 'sample_id', template=rtg_germline_vcf, axes_origin=[]),
            mgd.OutputFile('rtg_germline_maf', 'sample_id', template=rtg_germline_maf, axes_origin=[]),
            mgd.OutputFile('consensus_germline_maf', 'sample_id', template=consensus_germline_maf, axes_origin=[]),
            args['refdir'],
            normal_ids
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
