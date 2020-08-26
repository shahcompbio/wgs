import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import postprocessing
from wgs.workflows.postprocessing import tasks

def postprocessing_workflow(args):

    inputs = helpers.load_yaml(args['input_yaml'])

    samples = list(inputs.keys())

    #inputs
    pipeline_metadata = helpers.load_yaml(args['qc_metadata'])

    normal_bam_files = {sample: inputs[sample]['normal_bam'] for sample in samples}
    tumour_bam_files = {sample: inputs[sample]['tumour_bam'] for sample in samples}

    titan_files = {sample: inputs[sample]['titan'] for sample in samples}
    remixt_files = {sample: inputs[sample]['remixt'] for sample in samples}
    breakpoints_consensus_files = {sample: inputs[sample]['breakpoints_consensus'] for sample in samples}
    roh_files = {sample: inputs[sample]['roh'] for sample in samples}
    germline_calls_files = {sample: inputs[sample]['germline_calls'] for sample in samples}
    somatic_calls_files = {sample: inputs[sample]['somatic_calls'] for sample in samples}

    #outputs
    out_dir = args['out_dir']

    circos_plot_remixt = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_remixt.pdf')
    circos_plot_titan = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_titan.pdf')
    genome_wide_plot = os.path.join(out_dir, '{sample_id}', '{sample_id}_genome_wide.pdf')
    tumour_coverage = os.path.join(out_dir, '{sample_id}', '{sample_id}_tumour_coverage.tsv')
    normal_coverage = os.path.join(out_dir, '{sample_id}', '{sample_id}_normal_coverage.tsv')

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples
    )

    workflow.subworkflow(
        name="postprocessing",
        func=postprocessing.create_postprocessing_workflow,
        ctx=helpers.get_default_ctx(),
        axes=('sample_id',),
        args=(
            mgd.InputInstance('sample_id'),
            args["refdir"],
            pipeline_metadata,
            mgd.InputFile('normal.bam', 'sample_id', fnames=normal_bam_files),
            mgd.InputFile('tumour.bam', 'sample_id', fnames=tumour_bam_files),
            mgd.InputFile('titan', 'sample_id', fnames=titan_files),
            mgd.InputFile('remixt', 'sample_id', fnames=remixt_files),
            mgd.InputFile('breakpoints_consensus', 'sample_id', fnames=breakpoints_consensus_files),
            mgd.InputFile('roh', 'sample_id', fnames=roh_files),
            mgd.InputFile('germline_calls', 'sample_id', fnames=germline_calls_files),
            mgd.InputFile('somatic_calls', 'sample_id', fnames=somatic_calls_files),
            mgd.OutputFile('circos_plot_remixt.pdf', 'sample_id', template=circos_plot_remixt),
            mgd.OutputFile('circos_plot_titan.pdf', 'sample_id', template=circos_plot_titan),
            mgd.OutputFile('genome_wide_plot.pdf', 'sample_id', template=genome_wide_plot),
            mgd.OutputFile('normcov', 'sample_id', template=normal_coverage),
            mgd.OutputFile('tumcov', 'sample_id', template=tumour_coverage),
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list([circos_plot_remixt, circos_plot_titan, 
        genome_wide_plot], samples, "sample_id")
    meta_yaml = os.path.join(out_dir, 'pipeline_metadata.yaml')
    input_yaml_blob = os.path.join(out_dir, 'input.yaml')

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
            'metadata': {'type': 'postprocessing'}
        }
    )

    pyp.run(workflow)
