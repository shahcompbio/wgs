import os
import sys
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def get_output_paths(rootdir):
    return {"report": os.path.join(rootdir, "report.html"),
        "cna_table": os.path.join(rootdir, "cna_table.tsv"),
        "segments": os.path.join(rootdir, "segmental_copynumber.tsv"),
        "cohort_maf_oncogenic_filtered": os.path.join(rootdir, "cohort_oncogenic_filtered.maf")}


def load_inputs(yaml):
    yaml = helpers.load_yaml(yaml)
    cohort = list(yaml.keys())[0]
    data = yaml[cohort]
    outputs = {}
    for sample, sampledata in data.items():
        outputs[sample]  = sampledata
    return cohort, outputs


def cohort_qc_workflow(args):
    pypeline = pypeliner.app.Pypeline(config=args)

    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config.containers('wgs')))

    meta_yaml = os.path.join(args['out_dir'], 'metadata.yaml')
    input_yaml_blob = os.path.join(args['out_dir'], 'input.yaml')

    #inputs
    cohort, data = load_inputs(args['input_yaml'])

    out_dir = args["out_dir"]
    tmp_dir = args["tmpdir"]
    api_key = args["API_key"]
    metadata = helpers.load_yaml(os.path.join(args["refdir"], "metadata.yaml"))
    gtf = os.path.join(args["refdir"],metadata["paths"]["gtf"])

    germline_mafs = {label: data["germline_maf"] for label, data in data.items()}
    somatic_mafs = {label: data["somatic_maf"] for label, data in data.items()}
    remixt_data = {label: data["remixt"] for label, data in data.items()}

    #outputs
    filepaths = get_output_paths(os.path.join(out_dir, cohort))

    report_path = filepaths["report"]
    cna_table = filepaths["cna_table"]
    segmental_copynumber = filepaths["segments"]
    cohort_maf_oncogenic_filtered = filepaths["cohort_maf_oncogenic_filtered"]


    workflow.setobj(
        obj=mgd.OutputChunks('sample_label'),
        value=list(data.keys()),
    )


    workflow.subworkflow(
        name="classifycopynumber",
        func="wgs.workflows.cohort_qc.cna_annotation_workflow",
        args=(
            mgd.InputFile('remixt_dict', 'sample_label', fnames=remixt_data, axes_origin=[]),
            mgd.TempOutputFile('cna_maftools_table'),
            mgd.OutputFile(segmental_copynumber),
            mgd.OutputFile(cna_table),
            gtf,
        ),
    )

    workflow.subworkflow(
       name="maf_annotation_workflow",
       func="wgs.workflows.cohort_qc.preprocess_mafs_workflow",
       args=(
           mgd.InputFile('germline_mafs_dict', 'sample_label', fnames=germline_mafs, axes_origin=[]),
           mgd.InputFile('somatic_mafs_dict', 'sample_label', fnames=somatic_mafs, axes_origin=[]),
           mgd.OutputFile(cohort_maf_oncogenic_filtered),
           api_key
       ),
   )

    workflow.subworkflow(
       name="make_plots_and_report",
       func="wgs.workflows.cohort_qc.create_cohort_qc_report",
       args=(
           cohort,
           out_dir,
           mgd.InputFile(cohort_maf_oncogenic_filtered),
           mgd.TempInputFile('cna_maftools_table'),
           mgd.OutputFile(report_path)
       ),
   )   

    workflow.transform(
        name='generate_meta_files_results',
        func="wgs.utils.helpers.generate_and_upload_metadata",
        args=(
            sys.argv[0:],
            args["out_dir"],
            list(filepaths.values()),
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'cohort_qc'}
        }
    )

    pypeline.run(workflow)
