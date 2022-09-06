import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import breakpoint_calling_consensus
from wgs.workflows import destruct_wgs
from wgs.workflows import lumpy
from wgs.workflows import svaba


def breakpoint_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args["out_dir"], 'metadata.yaml')
    input_yaml_blob = os.path.join(args["out_dir"], 'input.yaml')

    tumour = inputs['tumour']
    normal = inputs['normal']
    sample_id = inputs['sample_id']

    destruct_breakpoints = args['output_prefix'] + '_destruct_breakpoints.csv.gz'
    destruct_library = args['output_prefix'] + '_destruct_library.csv.gz'
    destruct_raw_breakpoints = args['output_prefix'] + '_destruct_raw_breakpoints.csv.gz'
    destruct_raw_library = args['output_prefix'] + '_destruct_raw_library.csv.gz'
    destruct_reads = args['output_prefix'] + '_destruct_reads.csv.gz'

    lumpy_vcf = args['output_prefix'] + '_lumpy.vcf'

    svaba_vcf = args['output_prefix'] + '_svaba.vcf'

    parsed_csv = args['output_prefix'] + '_filtered_consensus_calls.csv.gz'

    single_node = args['single_node']

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='destruct',
        func=destruct_wgs.create_destruct_wgs_workflow,
        args=(
            mgd.InputFile(tumour, extensions=['.bai']),
            mgd.InputFile(normal, extensions=['.bai']),
            mgd.OutputFile(destruct_raw_breakpoints),
            mgd.OutputFile(destruct_raw_library),
            mgd.OutputFile(destruct_breakpoints),
            mgd.OutputFile(destruct_library),
            mgd.OutputFile(destruct_reads),
            sample_id,
            refdir_paths['reference'],
            refdir_paths['refdata_destruct'],
            refdir_paths['gtf']
        ),
        kwargs={'single_node': single_node}
    )

    workflow.subworkflow(
        name='lumpy',
        func=lumpy.create_lumpy_workflow,
        args=(
            mgd.OutputFile(lumpy_vcf),
        ),
        kwargs={
            'tumour_bam': mgd.InputFile(tumour, extensions=['.bai']),
            'normal_bam': mgd.InputFile(normal, extensions=['.bai']),
            'single_node': single_node
        },
    )

    if args['svaba']:
        workflow.subworkflow(
            name='svaba',
            func=svaba.create_svaba_workflow,
            args=(
                mgd.InputFile(tumour, extensions=['.bai']),
                mgd.InputFile(normal, extensions=['.bai']),
                mgd.OutputFile(svaba_vcf),
                refdir_paths['reference'],
            ),
        )

    workflow.subworkflow(
        name="consensus_calling",
        func=breakpoint_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile(destruct_breakpoints),
            mgd.InputFile(lumpy_vcf),
            mgd.OutputFile(parsed_csv, extensions=['.yaml']),
            chromosomes
        ),
    )

    filenames = [
        destruct_breakpoints,
        destruct_library,
        destruct_raw_breakpoints,
        destruct_raw_library,
        destruct_reads,
        lumpy_vcf,
        parsed_csv
    ]

    if args['svaba']:
        filenames.append(svaba_vcf)

    workflow.transform(
        name='generate_meta_files_results',
        func=helpers.generate_and_upload_metadata,
        args=(
            sys.argv[0:],
            args["out_dir"],
            filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'breakpoint_calling'}
        }
    )

    pyp.run(workflow)
