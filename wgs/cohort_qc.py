import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def cohort_qc_workflow(args):
    pypeline = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    inputs = helpers.load_qc_input_yaml_flat(args['input_yaml'])
    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]
    api_key = args["API_key"]

    sample_mafs = {label: data["sample_maf"] for label, data in inputs.items()}
    sample_labels = {label: data["sample_label"] for label, data in inputs.items()}

    report_path = {label[0]: os.path.join(out_dir, label[0], "report.html") for label, data in inputs.items()}

    workflow.setobj(
        obj=mgd.OutputChunks('cohort_label', 'sample_label'),
        value=list(inputs.keys()),
    )

    workflow.subworkflow(
        name="run_cohort_qc_workflow",
        func="wgs.workflows.cohort_qc.create_cohort_qc_workflow",
        axes=("cohort_label",),
        args=(
            mgd.InputInstance("cohort_label", ),
            api_key,
            out_dir,
            sample_labels, 
            mgd.InputFile('sample_maf_dict', 'cohort_label', 'sample_label', fnames=sample_mafs, axes_origin=[]),
            mgd.OutputFile('report_path', 'cohort_label', fnames=report_path),
        
        ),
    )

    pypeline.run(workflow)
