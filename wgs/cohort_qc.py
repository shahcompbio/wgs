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
    metadata = helpers.load_yaml(os.path.join(args["refdir"], "metadata.yaml"))
    gtf = os.path.join(args["refdir"],metadata["paths"]["gtf"])

    germline_mafs = {label: data["germline_maf"] for label, data in inputs.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in inputs.items()}
    remixt_data = {label: data["remixt"] for label, data in inputs.items()}
    sample_labels = {label: data["group_label"] for label, data in inputs.items()}

    report_path = {label[0]: os.path.join(out_dir, label[0], "report.html") for label, data in inputs.items()}
    cna_table = {label[0]: os.path.join(out_dir, label[0], "cna_table.tsv") for label, data in inputs.items()}
    segmental_copynumber = {label[0]: os.path.join(out_dir, label[0], "segmental_copynumber.tsv") for label, data in inputs.items()}
    cohort_maf = {label[0]: os.path.join(out_dir, label[0], "cohort.maf") for label, data in inputs.items()}
    cohort_maf_oncogenic_filtered = {label[0]: os.path.join(out_dir, label[0], "cohort_oncogenic_filtered.maf") for label, data in inputs.items()}
    cohort_cna = {label[0]: os.path.join(out_dir, label[0], "cohort_cna.tsv") for label, data in inputs.items()}


    workflow.setobj(
        obj=mgd.OutputChunks('cohort_label', 'sample_label'),
        value=list(inputs.keys()),
    )

    workflow.subworkflow(
        name="classifycopynumber",
        func="wgs.workflows.cohort_qc.cna_annotation_workflow",
        axes=("cohort_label",),
        args=(
            sample_labels,
            mgd.InputInstance('cohort_label'),
            mgd.InputFile('remixt_dict', 'cohort_label', 'sample_label', fnames=remixt_data, axes_origin=[]),
            mgd.OutputFile('cna_table', 'cohort_label', fnames=cna_table),
            mgd.OutputFile('segmental_copynumber', 'cohort_label', fnames=segmental_copynumber),
            gtf,
        ),
    )

    workflow.subworkflow(
        name="maf_annotation_workflow",
        func="wgs.workflows.cohort_qc.preprocess_mafs_workflow",
        axes=("cohort_label",),
        args=(
            mgd.InputInstance("cohort_label",),
            sample_labels,
            mgd.InputFile('germline_mafs_dict', 'cohort_label', 'sample_label', fnames=germline_mafs, axes_origin=[]),
            mgd.InputFile('somatic_mafs_dict', 'cohort_label', 'sample_label', fnames=somatic_mafs, axes_origin=[]),
            mgd.OutputFile('cohort_maf_oncogenic_filtered', 'cohort_label', fnames=cohort_maf_oncogenic_filtered),
            api_key
        ),
    )



    workflow.subworkflow(
        name="make_plots_and_report",
        func="wgs.workflows.cohort_qc.create_cohort_qc_report",
        axes=("cohort_label",),
        args=(
            mgd.InputInstance("cohort_label",),
            out_dir,
            mgd.InputFile('cohort_maf_oncogenic_filtered', 'cohort_label', fnames=cohort_maf_oncogenic_filtered),
            mgd.InputFile('cna_table', 'cohort_label', fnames=cna_table),
            mgd.OutputFile('report_path', 'cohort_label', fnames=report_path)
        ),
    )   


    pypeline.run(workflow)
