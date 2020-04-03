import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows import postprocessing


def postprocessing_workflow(args):
    yamldata = yaml.safe_load(open(args['input_yaml']))

    samples = yamldata.keys()

    normals = {sample: yamldata[sample]['normal'] for sample in samples}
    tumours = {sample: yamldata[sample]['tumour'] for sample in samples}

    variant_dir = {sample: yamldata[sample]['variant_dir'] for sample in samples}
    breakpoint_dir = {sample: yamldata[sample]['breakpoint_dir'] for sample in samples}
    copynumber_dir = {sample: yamldata[sample]['copynumber_dir'] for sample in samples}

    out_dir = args['out_dir']
    meta_yaml = os.path.join(out_dir, 'metadata.yaml')
    input_yaml_blob = os.path.join(out_dir, 'input.yaml')

    circos_plot = os.path.join(out_dir, '{sample_id}', '{sample_id}_circos.pdf')
    genome_wide_plot = os.path.join(out_dir, '{sample_id}', '{sample_id}_genome_wide.pdf')

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
            mgd.InputFile('normal.bam', 'sample_id', fnames=normals),
            mgd.InputFile('tumour.bam', 'sample_id', fnames=tumours),
            variant_dir,
            copynumber_dir,
            breakpoint_dir,
            mgd.OutputFile('circos_plot.pdf', 'sample_id', template=circos_plot),
            mgd.OutputFile('genome_wide_plot.pdf', 'sample_id', template=genome_wide_plot),
            args['refdir'],
            mgd.InputInstance('sample_id'),
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list([circos_plot, genome_wide_plot], samples, "sample_id")

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
            'metadata': {'type': 'copynumber_calling'}
        }
    )

    pyp.run(workflow)
