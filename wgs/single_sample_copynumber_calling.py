import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import hmmcopy


def single_sample_copynumber_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    inputs = helpers.load_yaml(args['input_yaml'])

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    bam = inputs['bam']
    sample_id = inputs['sample_id']

    bias_pdf = args['output_prefix'] + '_bias.pdf'
    correction_pdf = args['output_prefix'] + '_correction.pdf'
    hmmcopy_pdf = args['output_prefix'] + '_hmmcopy.pdf'
    correction_table = args['output_prefix'] + '_correctreads_with_state.txt'
    pygenes = args['output_prefix'] + '_hmmcopy.seg.pygenes'

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='hmmcopy',
        func=hmmcopy.create_hmmcopy_workflow,
        args=(
            mgd.InputFile(bam, extensions=['.bai']),
            sample_id,
            mgd.OutputFile(bias_pdf),
            mgd.OutputFile(correction_pdf),
            mgd.OutputFile(hmmcopy_pdf),
            mgd.OutputFile(correction_table),
            mgd.OutputFile(pygenes),
            chromosomes,
            refdir_paths['map_wig'],
            refdir_paths['gc_wig'],
            refdir_paths['gtf']
        ),
    )

    filenames = [
        bias_pdf,
        correction_pdf,
        hmmcopy_pdf,
        correction_table,
        pygenes,
    ]

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'single_sample_copynumber_calling'}
        }
    )

    pyp.run(workflow)
