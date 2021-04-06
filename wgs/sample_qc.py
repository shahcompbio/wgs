import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import sample_qc

def make_inputs(inputs, normal_only=False):
    samples = list(inputs.keys())
    roh = {sample: inputs[sample]['roh'] for sample in samples}
    germline = {sample: inputs[sample]['germline_calls'] for sample in samples}
    normal = {sample: inputs[sample]['normal_bam'] for sample in samples}

    tumor=None
    titan=None
    remixt=None
    breakpoints=None
    somatic=None

    if not normal_only:
        tumor = {sample: inputs[sample]['tumour_bam'] for sample in samples}
        titan = {sample: inputs[sample]['titan'] for sample in samples}
        remixt = {sample: inputs[sample]['remixt'] for sample in samples}
        breakpoints = {sample: inputs[sample]['breakpoints_consensus'] for sample in samples}
        somatic = {sample: inputs[sample]['somatic_calls'] for sample in samples}

    return {"normal":normal, "tumor": tumor, "roh":roh, "germline":germline, 
        "titan": titan, "remixt": remixt, "breakpoints": breakpoints, "somatic": somatic}


def sample_qc_workflow(args):
    inputs = helpers.load_yaml(args['input_yaml'])
    normal_only = args['normal_only']
    samples = list(inputs.keys())

    # inputs
    chromosomes = config.refdir_data(args['refdir'])['params']['chromosomes']
    files = make_inputs(inputs, normal_only=normal_only)

    # outputs
    out_dir = args['out_dir']
    normal_coverage = os.path.join(out_dir, '{sample_id}', '{sample_id}_normal_coverage.tsv')
    genome_wide_plot = os.path.join(out_dir, '{sample_id}', '{sample_id}_genome_wide.pdf')

    if not normal_only:
        circos_plot_remixt = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_remixt.pdf')
        circos_plot_titan = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos_titan.pdf')
        tumour_coverage = os.path.join(out_dir, '{sample_id}', '{sample_id}_tumour_coverage.tsv')


    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples
    )

    if normal_only:
        workflow.subworkflow(
            name="normal_sample_qc",
            func=sample_qc.create_sample_qc_workflow_normal_only,
            ctx=helpers.get_default_ctx(),
            axes=('sample_id',),
            args=(
                mgd.InputInstance('sample_id'),
                args["refdir"],
                mgd.InputFile('normal.bam', 'sample_id', fnames=files["normal"]),
                mgd.InputFile('roh', 'sample_id', fnames=files["roh"]),
                mgd.InputFile('germline_calls', 'sample_id', fnames=files["germline"]),
                mgd.OutputFile('genome_wide_plot.pdf', 'sample_id', template=genome_wide_plot),
                mgd.OutputFile('normcov', 'sample_id', template=normal_coverage),
                chromosomes,
                args['bins'],
                args['mapping_qual_threshold']
            ),
            # kwargs={'single_node': args['single_node']}
        )
        outputted_filenames = helpers.expand_list(
            [normal_coverage,genome_wide_plot],
            samples, "sample_id")
    else:
        workflow.subworkflow(
            name="sample_qc",
            func=sample_qc.create_sample_qc_workflow,
            ctx=helpers.get_default_ctx(),
            axes=('sample_id',),
            args=(
                mgd.InputInstance('sample_id'),
                args["refdir"],
                mgd.InputFile('normal.bam', 'sample_id', fnames=files["normal"]),
                mgd.InputFile('tumour.bam', 'sample_id', fnames=files["tumor"]),
                mgd.InputFile('titan', 'sample_id', fnames=files["titan"]),
                mgd.InputFile('remixt', 'sample_id', fnames=files["remixt"]),
                mgd.InputFile('breakpoints_consensus', 'sample_id', fnames=files["breakpoints"]),
                mgd.InputFile('roh', 'sample_id', fnames=files["roh"]),
                mgd.InputFile('germline_calls', 'sample_id', fnames=files["germline"]),
                mgd.InputFile('somatic_calls', 'sample_id', fnames=files["somatic"]),
                mgd.OutputFile('genome_wide_plot.pdf', 'sample_id', template=genome_wide_plot),
                mgd.OutputFile('normcov', 'sample_id', template=normal_coverage),
                mgd.OutputFile('tumcov', 'sample_id', template=tumour_coverage),
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
            axes=('sample_id',),
            func=sample_qc.circos_plot,
            args=(
                mgd.InputFile('titan', 'sample_id', fnames=files["titan"]),
                mgd.InputFile('remixt', 'sample_id', fnames=files["remixt"]),
                mgd.InputInstance("sample_id"),
                mgd.InputFile(
                    'breakpoints_consensus', 'sample_id', 
                    fnames=files["breakpoints"]),
                mgd.OutputFile(
                    'circos_remixt', 'sample_id', 
                    template=circos_plot_remixt),
                mgd.OutputFile('circos_titan', 'sample_id', template=circos_plot_titan),
            ),
        )
        outputted_filenames = helpers.expand_list(
            [circos_plot_remixt, circos_plot_titan, 
            normal_coverage, tumour_coverage,genome_wide_plot],
            samples, "sample_id")

    meta_yaml = os.path.join(out_dir, 'metadata.yaml')
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
            'metadata': {'type': 'sample_qc'}
        }
    )

    pyp.run(workflow)
