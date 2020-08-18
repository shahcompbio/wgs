import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
import os
import logging

def create_cohort_qc_workflow(cohort_label, api_key, out_dir, tmp_dir, cohort_maf,
    report_path
):

    #outputs
    oncoplot = os.path.join(out_dir, cohort_label, "cohort_oncoplot.png") 
    somatic_interactions_plot = os.path.join(out_dir, cohort_label, "somatic_interactions.png") 
    summary_plot = os.path.join(out_dir, cohort_label, "summary.png") 
    burden_plot = os.path.join(out_dir, cohort_label, "mutation_burden.png") 

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config.containers('wgs')})


    if api_key:
        filtered_maf = os.path.join(out_dir, cohort_label, "onco_kb-filtered_maf.maf") 
        annotated_maf = os.path.join(tmp_dir, cohort_label, "oncokb-annotated_maf.maf")

        workflow.transform(
            name='annotate_maf',
            func='wgs.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
            args=(
                cohort_maf,
                api_key,
                mgd.TempSpace("annotated_maf_tm"),
                mgd.TempOutputFile("annotated_maf"),
            ),
        )
        
        workflow.transform(
            name='filter_annotated_maf',
            func='wgs.workflows.cohort_qc.tasks.filter_annotated_maf',
            args=(
                mgd.TempInputFile("annotated_maf"),
                mgd.OutputFile(filtered_maf),
            ),
        )

        workflow.transform(
            name='make_cohort_plots',
            func='wgs.workflows.cohort_qc.tasks.make_R_cohort_plots',
            args=(
                cohort_maf,
                mgd.OutputFile(oncoplot),
                mgd.OutputFile(somatic_interactions_plot),
                mgd.OutputFile(summary_plot),
            ),
            kwargs={
                "filtered_maf":mgd.InputFile(filtered_maf)
            },
        )

    if not api_key:
        logging.warning("No API key is provided to use oncoKB, so results will be unfiltered.")

        workflow.transform(
            name='make_R_cohort_plots',
            func='wgs.workflows.cohort_qc.tasks.make_cohort_plots',
            args=(
                cohort_maf,
                mgd.OutputFile(oncoplot),
                mgd.OutputFile(somatic_interactions_plot),
                mgd.OutputFile(summary_plot),
            ),
        )


    workflow.transform(
        name='burden_plot',
        func='wgs.workflows.cohort_qc.tasks.plot_mutation_burden',
        args=(
            cohort_maf,
            mgd.OutputFile(burden_plot),
        ),
    )

    workflow.transform(
        name='make_report',
        func='wgs.workflows.cohort_qc.tasks.make_report',
        args=(
            cohort_label,
            mgd.InputFile(oncoplot),
            mgd.InputFile(somatic_interactions_plot),
            mgd.InputFile(summary_plot),
            mgd.InputFile(burden_plot),
            report_path,
        ),
    )

    return workflow
    