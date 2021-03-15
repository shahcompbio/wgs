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

    bams = helpers.get_values_from_input(inputs, 'bam')
    samples = list(bams.keys())

    cna_outdir = os.path.join(args['out_dir'], 'copynumber', '{sample_id}')

    hmmcopy_raw_dir = os.path.join(cna_outdir, 'hmmcopy')
    bias_pdf = os.path.join(hmmcopy_raw_dir, 'plots', '{sample_id}_bias.pdf')
    correction_pdf = os.path.join(hmmcopy_raw_dir, 'plots', '{sample_id}_correction.pdf')
    hmmcopy_pdf = os.path.join(hmmcopy_raw_dir, 'plots', '{sample_id}_hmmcopy.pdf')
    correction_table = os.path.join(hmmcopy_raw_dir, '{sample_id}_correctreads_with_state.txt')
    pygenes = os.path.join(hmmcopy_raw_dir, '{sample_id}_hmmcopy.seg.pygenes')

    refdir_paths = config.refdir_data(args['refdir'])['paths']
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='hmmcopy',
        func=hmmcopy.create_hmmcopy_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("sample.bam", 'sample_id', fnames=bams,
                          extensions=['.bai']),
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('bias', 'sample_id', template=bias_pdf),
            mgd.OutputFile('correction', 'sample_id', template=correction_pdf),
            mgd.OutputFile('hmmcopy', 'sample_id', template=hmmcopy_pdf),
            mgd.OutputFile('correction_table', 'sample_id', template=correction_table),
            mgd.OutputFile('pygenes', 'sample_id', template=pygenes),
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

    outputted_filenames = helpers.expand_list(filenames, samples, "sample_id")

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            outputted_filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'single_sample_copynumber_calling'}
        }
    )

    pyp.run(workflow)
