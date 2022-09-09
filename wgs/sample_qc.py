import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import sample_qc


def sample_qc_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])

    roh = inputs['roh']
    germline = inputs['germline_calls']
    normal = inputs['normal_bam']
    tumour = inputs.get('tumour_bam', None)
    titan = inputs.get('titan', None)
    remixt = inputs.get('remixt', None)
    breakpoints = inputs.get('breakpoints_consensus', None)
    somatic = inputs.get('somatic_calls', None)
    sample_id = inputs['sample_id']

    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']

    # outputs
    normal_coverage = args['output_prefix'] + '_normal_coverage.tsv'
    genome_wide_plot = args['output_prefix'] + '_genome_wide.pdf'

    circos_plot_remixt = args['output_prefix'] + '_circos_remixt.pdf'
    circos_plot_titan = args['output_prefix'] + '_circos_titan.pdf'
    tumour_coverage = args['output_prefix'] + '_tumour_coverage.tsv'

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    if args['normal_only']:
        workflow.subworkflow(
            name="normal_sample_qc",
            func=sample_qc.create_sample_qc_workflow_normal_only,
            ctx=helpers.get_default_ctx(),
            args=(
                sample_id,
                args["refdir"],
                mgd.InputFile(normal),
                mgd.InputFile(roh),
                mgd.InputFile(germline),
                mgd.OutputFile(genome_wide_plot),
                mgd.OutputFile(normal_coverage),
                chromosomes,
                args['bins'],
                args['mapping_qual_threshold']
            ),
        )
        outputted_filenames = [normal_coverage, genome_wide_plot]
    else:
        workflow.subworkflow(
            name="sample_qc",
            func=sample_qc.create_sample_qc_workflow,
            ctx=helpers.get_default_ctx(),
            args=(
                sample_id,
                args["refdir"],
                mgd.InputFile(normal),
                mgd.InputFile(tumour),
                mgd.InputFile(titan),
                mgd.InputFile(remixt),
                mgd.InputFile(breakpoints),
                mgd.InputFile(roh),
                mgd.InputFile(germline),
                mgd.InputFile(somatic),
                mgd.OutputFile(genome_wide_plot),
                mgd.OutputFile(normal_coverage),
                mgd.OutputFile(tumour_coverage),
                chromosomes,
                args['bins'],
                args['mapping_qual_threshold']
            ),
            kwargs={'single_node': args['single_node']}
        )
        workflow.subworkflow(
            name='generate_circos_plot',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='24:00',
                disk=400
            ),
            func=sample_qc.circos_plot,
            args=(
                mgd.InputFile(titan),
                mgd.InputFile(remixt),
                sample_id,
                mgd.InputFile(breakpoints),
                mgd.OutputFile(circos_plot_remixt),
                mgd.OutputFile(circos_plot_titan),
            ),
        )
        outputted_filenames = [
            circos_plot_remixt, circos_plot_titan,
            normal_coverage, tumour_coverage, genome_wide_plot
        ]

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

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
            'metadata': {'type': 'sample_qc'}
        }
    )

    pyp.run(workflow)
