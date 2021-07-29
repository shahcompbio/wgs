import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def mutect_workflow(args):
    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    tumour = args['tumour']
    normal = args['normal']
    tumour_id = args['tumour_id']
    normal_id = args['normal_id']

    refdir = args['refdir']
    single_node = args['single_node']

    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    paths_refdir = config.refdir_data(refdir)['paths']

    var_dir = os.path.join(args['out_dir'], 'mutect')

    mutect_vcf = os.path.join(var_dir, 'mutect.vcf.gz')
    mutect_maf = os.path.join(var_dir, 'mutect.maf')

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name="mutect",
        func='wgs.workflows.mutect.create_mutect_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile(normal),
            mgd.InputFile(tumour),
            mgd.OutputFile(mutect_vcf),
            mgd.OutputFile(mutect_maf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
            tumour_id,
            normal_id
        ),
        kwargs={
            'single_node': single_node,
        },
    )

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args['out_dir'],
            [mutect_vcf, mutect_maf],
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'metadata': {'type': 'variant_calling'}
        }
    )

    pyp.run(workflow)
