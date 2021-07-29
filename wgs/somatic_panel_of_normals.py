import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import somatic_panel_of_normals


def somatic_panel_of_normals_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    inputs = inputs['normals']
    samples = list(inputs.keys())
    tumours = {samp: inputs[samp] for samp in samples}

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    var_dir = os.path.join(args['out_dir'], 'somatic_panel_of_normals')
    pon_vcf = os.path.join(var_dir, 'mutect2_panel_of_normals.vcf.gz')
    sample_vcf = os.path.join(var_dir, '{sample_id}_mutect.vcf.gz')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    workflow.subworkflow(
        name='somatic_panel_of_normals',
        func=somatic_panel_of_normals.create_somatic_panel_of_normals_workflow,
        args=(
            list(samples),
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile(pon_vcf),
            mgd.OutputFile('mutect_vcf', 'sample_id', template=sample_vcf, axes_origin=[]),
            args['refdir'],
            args['single_node']
        ),
    )

    outputted_filenames = helpers.expand_list([pon_vcf, sample_vcf], samples, "sample_id")

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
