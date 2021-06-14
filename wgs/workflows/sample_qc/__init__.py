import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows.sample_qc import tasks

def circos_plot(titan_calls, remixt_calls, sample_id, breakpoints,
           circos_plot_remixt, circos_plot_titan):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='prep_titan',
        func='wgs_qc_utils.reader.read_titan.make_for_circos',
        ctx=helpers.get_default_ctx(
            memory=5
        ),
        args=(
            mgd.InputFile(titan_calls),
            mgd.TempOutputFile("titan_prepped"),
        )
    )

    workflow.transform(
        name='prep_remixt',
        func='wgs_qc_utils.reader.read_remixt.make_for_circos',
        ctx=helpers.get_default_ctx(
            memory=5
        ),
        args=(
            mgd.InputFile(remixt_calls),
            sample_id,
            mgd.TempOutputFile("remixt_prepped"),
        )
    )
    workflow.transform(
        name='circos_plot',
        func='wgs.workflows.sample_qc.tasks.circos',
        ctx=helpers.get_default_ctx(
            memory=5
        ),
        args=(
            mgd.TempInputFile("titan_prepped"),
            mgd.TempInputFile("remixt_prepped"),
            sample_id,
            breakpoints,
            mgd.OutputFile(circos_plot_remixt),
            mgd.OutputFile(circos_plot_titan),
            mgd.TempSpace("circos")
        )
    )

    return workflow

def get_coverage_data(
        input_bam, output, refdir, chromosomes,
        mapping_qual, bins, single_node=False
):
    reference = config.refdir_data(refdir)['paths']['reference']

    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='generate_coverage_bed',
            func='wgs.workflows.sample_qc.tasks.generate_coverage_bed',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                reference,
                mgd.TempOutputFile('coverage_bed.bed'),
                chromosomes,
                bins,
            )
        )
        workflow.transform(
            name='samtools_coverage',
            func='wgs.workflows.sample_qc.tasks.samtools_coverage',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                mapping_qual,

            ),
        )

    else:
        workflow.setobj(
            obj=mgd.OutputChunks('chromosome'),
            value=chromosomes
        )
        workflow.transform(
            name='generate_coverage_bed',
            func='wgs.workflows.sample_qc.tasks.generate_coverage_bed',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            axes=('chromosome',),
            args=(
                reference,
                mgd.TempOutputFile('coverage_bed.bed', 'chromosome'),
                mgd.InputInstance('chromosome'),
                bins,
            )
        )
        workflow.transform(
            name='samtools_coverage',
            func='wgs.workflows.sample_qc.tasks.samtools_coverage',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            axes=('chromosome',),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed', 'chromosome'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                mapping_qual,
            ),
        )

        workflow.transform(
            name='merge_data',
            func='wgs.utils.csvutils.concatenate_csv',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                mgd.TempInputFile('per_interval.txt', 'chromosome', axes_origin=[]),
                mgd.OutputFile(output),
            )
        )

    return workflow




def create_sample_qc_workflow(
        sample_id,
        refdir,
        normal_bam,
        tumour_bam,
        titan,
        remixt,
        breakpoints_consensus,
        roh,
        germline_calls,
        somatic_calls,
        genome_wide_plot,
        normal_coverage,
        tumour_coverage,
        chromosomes,
        bins,
        mapping_qual_threshold,
        single_node=False
):

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='coverage_normal_data',
        func=get_coverage_data,
        args=(
            mgd.InputFile(normal_bam),
            mgd.OutputFile(normal_coverage),
            refdir,
            chromosomes,
            mapping_qual_threshold,
            bins,
        ),
        kwargs={'single_node': single_node}
    )

    workflow.subworkflow(
        name='coverage_tumour_data',
        func=get_coverage_data,
        args=(
            mgd.InputFile(tumour_bam),
            mgd.OutputFile(tumour_coverage),
            refdir,
            chromosomes,
            mapping_qual_threshold,
            bins,
        ),
        kwargs={'single_node': single_node}
    )


    workflow.transform(
        name='generate_genome_wide_plot',
        ctx=helpers.get_default_ctx(
            memory=10,
        ),
        func="wgs.workflows.sample_qc.tasks.genome_wide",
        args=(
            sample_id,
            mgd.InputFile(roh),
            mgd.InputFile(germline_calls),
            mgd.InputFile(normal_coverage),
            chromosomes,
            mgd.OutputFile(genome_wide_plot),
        ),
        kwargs={"titan": mgd.InputFile(titan),
            "somatic": mgd.InputFile(somatic_calls),
            "remixt": mgd.InputFile(remixt),
            "tumour": mgd.InputFile(tumour_coverage),
            "breakpoints": mgd.InputFile(breakpoints_consensus)
        }
    )

    return workflow


def create_sample_qc_workflow_normal_only(
        sample_id,
        refdir,
        normal_bam,
        roh,
        germline_calls,
        genome_wide_plot,
        normal_coverage,
        chromosomes,
        bins,
        mapping_qual_threshold,
        single_node=False
):

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='coverage_normal_data',
        func=get_coverage_data,
        args=(
            mgd.InputFile(normal_bam),
            mgd.OutputFile(normal_coverage),
            refdir,
            chromosomes,
            mapping_qual_threshold,
            bins,
        ),
        kwargs={'single_node': single_node}
    )



    workflow.transform(
        name='generate_genome_wide_plot',
        ctx=helpers.get_default_ctx(
            memory=10,
        ),
        func="wgs.workflows.sample_qc.tasks.genome_wide",
        args=(
            sample_id,
            mgd.InputFile(roh),
            mgd.InputFile(germline_calls),
            mgd.InputFile(normal_coverage),
            chromosomes,
            mgd.OutputFile(genome_wide_plot),
        ),
        kwargs={"normal_only":True}
    )

    return workflow
