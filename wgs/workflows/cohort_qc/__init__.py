import logging
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


def cna_annotation_workflow(remixt_dict, output_table, segmental_copynumber, cbio_cna_table, gtf):
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
            mgd.TempOutputFile('cn_change', 'sample_label'),
        ),
    )

    workflow.transform(
        name='merge_cn_tables',
        func='wgs.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('cn_change', 'sample_label', axes_origin=[]),
            mgd.TempOutputFile("merged_cn_change"),
        ),
    )

    workflow.transform(
        name='make_cbio_cna_table',
        func='wgs.workflows.cohort_qc.tasks.make_cbio_cna_table',
        args=(
            mgd.TempInputFile('merged_cn_change'),
            mgd.OutputFile(cbio_cna_table),
        ),
    )

    workflow.transform(
        name='make_maftools_cna_table',
        func='wgs.workflows.cohort_qc.tasks.make_maftools_cna_table',
        args=(
            mgd.TempInputFile('merged_cn_change'),
            output_table
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


def preprocess_mafs_workflow(germline_maf_dict, somatic_maf_dict, merged_annotated_maf, api_key
):

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label'),
        value=list(germline_maf_dict.keys()),
    )

    workflow.transform(
        name='annotate_germline_mafs',
        func='wgs.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile('germlne_maf', 'sample_label', fnames=germline_maf_dict),
            api_key,
            mgd.TempSpace("annotated_germline_maf_tmp", 'sample_label'),
            mgd.TempOutputFile("annotated_germline_maf", 'sample_label'),
        ),
    )
    workflow.transform(
        name='filter_germline_mafs',
        func='wgs.workflows.cohort_qc.tasks.filter_maf',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile("annotated_germline_maf", 'sample_label'),
            mgd.TempOutputFile("filtered_germline_maf", "sample_label")
        ),
    )
    workflow.transform(
        name='annotate_somatic_mafs',
        func='wgs.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile('somatic_maf', 'sample_label', fnames=somatic_maf_dict),
            api_key,
            mgd.TempSpace("annotated_somatic_maf_tmp", 'sample_label'),
            mgd.TempOutputFile("annotated_somatic_maf", 'sample_label'),
        ),
    )

    workflow.transform(
        name='annotate_germline_class',
        func='wgs.workflows.cohort_qc.tasks.annotate_germline_somatic',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile('filtered_germline_maf', 'sample_label'),
            mgd.TempOutputFile("filtered_class_labeled_germline_maf", 'sample_label'),
            True
        ),
    )

    workflow.transform(
        name='annotate_somatic_class',
        func='wgs.workflows.cohort_qc.tasks.annotate_germline_somatic',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile('annotated_somatic_maf', 'sample_label'),
            mgd.TempOutputFile("class_labeled_somatic_maf", 'sample_label'),
            False
        ),
    )

    workflow.transform(
        name='merge_filtered_germline_somatic',
        func='wgs.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile('filtered_class_labeled_germline_maf', 'sample_label', axes_origin=[]),
            mgd.TempInputFile('class_labeled_somatic_maf', 'sample_label', axes_origin=[]),
            mgd.OutputFile(merged_annotated_maf),
        ),
    )

    return workflow


def create_cohort_qc_report(
    cohort_label, out_dir, filtered_cohort_maf, cna_table, report_path
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

    non_synonymous_labels=["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
        "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
        "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
    ]

    workflow.transform(
        name='postprocess_maf',
        func='wgs.workflows.cohort_qc.tasks.prepare_maf_for_maftools',
        args=(
            cohort_label,
            mgd.InputFile(filtered_cohort_maf),
            mgd.TempOutputFile("prepared_maf"),
            non_synonymous_labels,
            mgd.TempOutputFile("vcNames")
        ),
    )


    workflow.transform(
        name='burden_plot',
        func='wgs.workflows.cohort_qc.tasks.plot_mutation_burden',
        args=(
            mgd.InputFile(filtered_cohort_maf),
            mgd.OutputFile(burden_plot),
        ),
    )

    workflow.transform(
        name='build_gene_list',
        func='wgs.workflows.cohort_qc.tasks.build_gene_list',
        args=(
            mgd.InputFile(cna_table),
            mgd.TempOutputFile("genelist")
        ),
    )
    workflow.transform(
        name='make_cohort_plots',
        func='wgs.workflows.cohort_qc.tasks.make_R_cohort_plots',
        args=(
            mgd.TempInputFile("prepared_maf"),
            mgd.InputFile(cna_table),
            mgd.OutputFile(oncoplot),
            mgd.OutputFile(somatic_interactions_plot),
            mgd.OutputFile(summary_plot),
            mgd.TempInputFile("vcNames"),
            mgd.TempInputFile("genelist")

        ),
        kwargs={'docker_image':config.containers("wgs_qc_html") },
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
