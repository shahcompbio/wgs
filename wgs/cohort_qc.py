import os
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def get_cbioportal_paths(root_dir):
    """Get cbioportal output paths.

    Args:
        root_dir ([str]): [path to out_dir]

    Returns:
        [dict]: [labeled output paths]
    """
    if not os.path.exists(os.path.join(root_dir, "cbioportal")):
        os.makedirs(os.path.join(root_dir, "cbioportal"))

    filtered_germline_maf = os.path.join(
        root_dir, "cbioportal", "filtered_germline.maf"
    )
    annotated_somatic_maf = os.path.join(
        root_dir, "cbioportal", "annotated_somatic.maf"
    )
    cna_table = os.path.join(
        root_dir, "cbioportal",  "cna_table.tsv"
    )
    segments = os.path.join(
        root_dir, "cbioportal",  "segments.tsv"
    )

    return {
        "filtered_germline_maf": filtered_germline_maf,
        "annotated_somatic_maf": annotated_somatic_maf,
        "cna_table": cna_table,
        "segments": segments
    }


def get_maftools_paths(root_dir):
    """Get maftools output paths.

    Args:
        root_dir ([str]): [path to out_dir]

    Returns:
        [dict]: [labeled output paths]
    """

    if not os.path.exists(os.path.join(root_dir, "maftools")):
        os.makedirs(os.path.join(root_dir, "maftools"))

    cohort_oncoplot = os.path.join(
        root_dir, "maftools", "cohort_oncoplot.png"
    )
    maftools_maf = os.path.join(
        root_dir,  "maftools", "maftools_maf.maf"
    )
    maftools_cna = os.path.join(
        root_dir, "maftools",  "maftools_cna.tsv"
    )
    report = os.path.join(
        root_dir, "maftools",  "report.html"
    )

    return {
        "cohort_oncoplot": cohort_oncoplot,
        "maftools_maf": maftools_maf,
        "maftools_cna": maftools_cna,
        "report": report
    }




def cohort_qc_workflow(args):
    pypeline = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    cohort, data = helpers.load_cohort_qc_inputs(args['input_yaml'])

    out_dir = args["out_dir"]
    api_key = args["API_key"]
    metadata = helpers.load_yaml(os.path.join(args["refdir"], "metadata.yaml"))
    gtf = os.path.join(args["refdir"], metadata["paths"]["gtf"])

    germline_mafs = {label: data["germline_maf"] for label, data in data.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in data.items()}
    remixt_data = {label: data["remixt"] for label, data in data.items()}

    maftools_paths = get_maftools_paths(os.path.join(out_dir, cohort))
    cbio_paths = get_cbioportal_paths(os.path.join(out_dir, cohort))

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label'),
        value=list(data.keys()),
    )



    workflow.subworkflow(
        name="classifycopynumber",
        func="wgs.workflows.cohort_qc.cna_annotation_workflow",
        args=(
            mgd.InputFile(
                'remixt_dict', 'sample_label',
                fnames=remixt_data, axes_origin=[]
            ),
            mgd.OutputFile(cbio_paths["cna_table"]),
            mgd.OutputFile(maftools_paths["maftools_cna"]),
            mgd.OutputFile(cbio_paths["segments"]),
            gtf
        ),
    )


    workflow.subworkflow(
        name="maf_annotation_workflow",
        func="wgs.workflows.cohort_qc.preprocess_mafs_workflow",
        args=(
            mgd.InputFile(
                'germline_mafs_dict','sample_label',
                fnames=germline_mafs, axes_origin=[]
            ),
            mgd.InputFile(
                'somatic_mafs_dict', 'sample_label', 
                fnames=somatic_mafs, axes_origin=[]
            ),
            mgd.OutputFile(cbio_paths["filtered_germline_maf"]),
            mgd.OutputFile(cbio_paths["annotated_somatic_maf"]),
            api_key
        ),
    )



    workflow.subworkflow(
        name="make_plots_and_report",
        func="wgs.workflows.cohort_qc.create_cohort_qc_report",
        args=(
            cohort,
            out_dir,
            mgd.InputFile(cbio_paths["filtered_germline_maf"]),
            mgd.InputFile(cbio_paths["annotated_somatic_maf"]),
            mgd.InputFile(maftools_paths["maftools_cna"]),
            mgd.InputFile(maftools_paths["maftools_maf"]),
            mgd.OutputFile(maftools_paths["report"]),
            mgd.OutputFile(maftools_paths["cohort_oncoplot"]),
        ),
    )

    pypeline.run(workflow)
