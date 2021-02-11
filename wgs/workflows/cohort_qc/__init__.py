import os
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def cna_annotation_workflow(
    remixt_dict,
    cbio_cna,
    maftools_cna,
    segments,
    gtf
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_label'),
        value=list(remixt_dict.keys()),
    )

    workflow.transform(
        name='classify_remixt',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
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
        name='generate_segmental_copynumber',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
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
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.merge_segmental_cn',
        args=(
            mgd.TempInputFile('segmental_cn', 'sample_label', axes_origin=[]),
            mgd.OutputFile(segments)
        ),
    )

    workflow.transform(
        name='merge_amp_tables',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('amps', 'sample_label', axes_origin=[]),
            mgd.TempOutputFile("merged_amps"),
        ),
    )

    workflow.transform(
        name='merge_del_tables',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.merge_cna_tables',
        args=(
            mgd.TempInputFile('dels', 'sample_label', axes_origin=[]),
            mgd.TempOutputFile("merged_dels"),
        ),
    )

    workflow.transform(
        name='make_cbio_cna_table',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.make_cbio_cna_table',
        args=(
            mgd.TempInputFile('merged_amps'),
            mgd.TempInputFile('merged_dels'),
            mgd.OutputFile(cbio_cna)
        ),
    )

    workflow.transform(
        name='make_maftools_cna_table',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.make_maftools_cna_table',
        args=(
            mgd.TempInputFile('merged_amps'),
            mgd.TempInputFile('merged_dels'),
            mgd.OutputFile(maftools_cna)
        ),
    )

    return workflow


def preprocess_mafs_workflow(
        germline_maf_dict,
        somatic_maf_dict,
        merged_germline,
        merged_somatic,
        api_key
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
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.annotate_maf_with_oncokb',
        axes=("sample_label",),
        args=(
            mgd.InputFile(
                'germlne_maf', 'sample_label', fnames=germline_maf_dict
            ),
            api_key,
            mgd.TempSpace("annotated_germline_maf_tmp", 'sample_label'),
            mgd.TempOutputFile("annotated_germline_maf", 'sample_label'),
        ),
    )

    workflow.transform(
        name='filter_germline_mafs',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.filter_maf',
        axes=("sample_label",),
        args=(
            mgd.TempInputFile("annotated_germline_maf", 'sample_label'),
            mgd.TempOutputFile("filtered_germline_maf", "sample_label")
        ),
    )

    workflow.transform(
        name='annotate_somatic_mafs',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
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
        name='merge_germline_mafs',
        func='wgs.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile(
                "filtered_germline_maf", 'sample_label',  axes_origin=[]
            ),
            mgd.OutputFile(merged_germline)
        ),
    )
    workflow.transform(
        name='merge_somatic_mafs',
        func='wgs.workflows.cohort_qc.tasks.merge_mafs',
        args=(
            mgd.TempInputFile(
                "annotated_somatic_maf", 'sample_label', axes_origin=[]
            ),
            mgd.OutputFile(merged_somatic)
        ),
    )

    return workflow



def create_cohort_qc_report(
        cohort_label,
        out_dir,
        filtered_germline,
        annotated_somatic,
        maftools_cna,
        maftools_maf,
        report_path,
        cohort_oncoplot
):
    somatic_interactions_plot = os.path.join(
        out_dir, cohort_label, "somatic_interactions.png"
    )
    summary_plot = os.path.join(out_dir, cohort_label, "summary.png")
    burden_plot = os.path.join(out_dir, cohort_label, "mutation_burden.png")

    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    ## TODO: move these to config
    non_synonymous_labels = [
        "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
        "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
        "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
    ]

    workflow.transform(
        name='postprocess_maf',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.prepare_maf_for_maftools',
        args=(
            mgd.InputFile(filtered_germline),
            mgd.InputFile(annotated_somatic),
            mgd.OutputFile(maftools_maf),
            non_synonymous_labels,
            mgd.TempOutputFile("vcNames")
        ),
    )


    workflow.transform(
        name='burden_plot',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.plot_mutation_burden',
        args=(
            mgd.InputFile(maftools_maf),
            mgd.OutputFile(burden_plot),
        ),
    )

    workflow.transform(
        name='build_gene_list',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='6:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.build_gene_list',
        args=(
            mgd.InputFile(maftools_cna),
            mgd.TempOutputFile("genelist")
        ),
    )

    workflow.transform(
        name='make_cohort_plots',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.make_R_cohort_plots',
        args=(
            mgd.InputFile(maftools_maf),
            mgd.InputFile(maftools_cna),
            mgd.OutputFile(cohort_oncoplot),
            mgd.OutputFile(somatic_interactions_plot),
            mgd.OutputFile(summary_plot),
            mgd.TempInputFile("vcNames"),
            mgd.TempInputFile("genelist")

        ),
        kwargs={'docker_image': config.containers("wgs_qc_html")},
    )

    workflow.transform(
        name='make_report',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func='wgs.workflows.cohort_qc.tasks.make_report',
        args=(
            cohort_label,
            mgd.InputFile(cohort_oncoplot),
            mgd.InputFile(somatic_interactions_plot),
            mgd.InputFile(summary_plot),
            mgd.InputFile(burden_plot),
            mgd.OutputFile(report_path),
        ),
        kwargs={'docker_image': config.containers("wgs_qc_html")}
    )

    return workflow
