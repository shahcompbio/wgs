import logging
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


def cna_annotation_workflow(cohort, labels, remixt_dict, output_table, segmental_copynumber, gtf):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )
    workflow.setobj(
        obj=mgd.OutputChunks('sample_label'),
        value=list(remixt_dict.keys()),
    )

    workflow.transform(
        name='classify_remixt',
        func='wgs.workflows.cohort_qc.tasks.classify_remixt',
        axes=("sample_label",),
        args=(
            mgd.InputInstance('sample_label'),
            mgd.InputFile('remixt', 'sample_label', fnames=remixt_dict),
            gtf,
            mgd.TempSpace('annotated_maf_tmp', 'sample_label'),
            mgd.TempOutputFile('amps', 'sample_label'),
            mgd.TempOutputFile('dels', 'sample_label'),
        ),
    )

    workflow.transform(
        name='merge_cna_tables',
        func='wgs.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('amps', 'sample_label', axes_origin=[]),
            mgd.TempInputFile('dels', 'sample_label', axes_origin=[]),
            labels,
            mgd.OutputFile(output_table),
            cohort
        ),
    )

    workflow.transform(
        name='generate_segmental_copynumber',
        func='wgs.workflows.cohort_qc.tasks.generate_segmental_copynumber',
        axes=("sample_label",),
        args=(
            mgd.InputFile('remixt', 'sample_label', fnames=remixt_dict),
            mgd.TempOutputFile('segmental_cn', 'sample_label'),
            mgd.InputInstance('sample_label')
        ),
    )

    workflow.transform(
        name='merge_segmental_cn',
        func='wgs.workflows.cohort_qc.tasks.merge_segmental_cn',
        args=(
            mgd.TempInputFile('segmental_cn', 'sample_label', axes_origin=[]),
            mgd.OutputFile(segmental_copynumber)
        ),
    )

    return workflow

def create_cohort_qc_workflow(
        cohort_label, api_key, out_dir, sample_labels, germline_mafs, somatic_mafs, 
        cna_table, report_path, cohort_maf_oncogenic_annotated, cohort_remixt_path
):

    oncoplot = os.path.join(
        out_dir, cohort_label, "cohort_oncoplot.png"
    )
    somatic_interactions_plot = os.path.join(
        out_dir, cohort_label, "somatic_interactions.png"
    )
    summary_plot = os.path.join(
        out_dir, cohort_label, "summary.png"
    )
    burden_plot = os.path.join(
        out_dir, cohort_label, "mutation_burden.png"
    )

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )


    workflow.transform(
        name='merge_germline_mafs',
        func='wgs.workflows.cohort_qc.tasks.merge_relabel_mafs',
        args=(
            germline_mafs,
            mgd.TempOutputFile("cohort_germline_maf"),
            sample_labels
        ),
        kwargs={"class_label": "_germline"}

    )

    workflow.transform(
        name='merge_somatic_mafs',
        func='wgs.workflows.cohort_qc.tasks.merge_relabel_mafs',
        args=(
            somatic_mafs,
            mgd.TempOutputFile("cohort_somatic_maf"),
            sample_labels
        ),
        kwargs={ "class_label": "_somatic"}

    )

    workflow.transform(
        name='merge_mafs',
        func='wgs.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            [mgd.TempInputFile("cohort_somatic_maf"),
            mgd.TempInputFile("cohort_germline_maf")],
            mgd.TempOutputFile("cohort_maf"),
        ),

    )
    
    if api_key:
        workflow.transform(
            name='annotate_maf',
            func='wgs.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
            args=(
                mgd.TempInputFile("cohort_maf"),
                api_key,
                mgd.TempSpace("annotated_maf_tmp"),
                mgd.OutputFile(cohort_maf_oncogenic_annotated),
            ),
            kwargs={'docker_image': config.containers("oncokb-annotator")}
        )

        workflow.transform(
            name='filter_annotated_maf',
            func='wgs.workflows.cohort_qc.tasks.filter_annotated_maf',
            args=(
                mgd.InputFile(cohort_maf_oncogenic_annotated),
                mgd.TempOutputFile("filtered_maf"),
                mgd.TempOutputFile('vcNames'),
            ),
        )
        kwargs={"filtered_maf":  mgd.TempInputFile("filtered_maf"), 'vcNames':mgd.TempInputFile('vcNames'),
            'docker_image':config.containers("wgs_qc_html") 
        }
    else:
        logging.warning("No API key is provided to use oncoKB, so results will be unfiltered.")
        kwargs = {'docker_image':config.containers("wgs_qc_html") }


    workflow.transform(
        name='make_cohort_plots',
        func='wgs.workflows.cohort_qc.tasks.make_R_cohort_plots',
        args=(
            mgd.InputFile(cohort_maf_oncogenic_annotated),
            mgd.InputFile(cna_table),
            mgd.OutputFile(oncoplot),
            mgd.OutputFile(somatic_interactions_plot),
            mgd.OutputFile(summary_plot),
        ),
        kwargs=kwargs,
    )

    workflow.transform(
        name='burden_plot',
        func='wgs.workflows.cohort_qc.tasks.plot_mutation_burden',
        args=(
            mgd.TempInputFile("cohort_somatic_maf"),
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
            mgd.OutputFile(report_path),
        ),
        kwargs={'docker_image':config.containers("wgs_qc_html") }
    )

    return workflow
