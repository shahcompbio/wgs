import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
import os

def get_coverage_data(input_bam, output, refdir, single_node=False):
    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    bins = config.refdir_data(refdir)['params']['bins']

    reference = config.refdir_data(refdir)['paths']['reference']

    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='generate_coverage_bed',
            func='wgs.workflows.postprocessing.tasks.generate_coverage_bed',
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
            func='wgs.workflows.postprocessing.tasks.samtools_coverage',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),

            ),
            kwargs={'docker_image': config.containers('samtools')},
        )

    else:

        workflow.setobj(
            obj=mgd.OutputChunks('chromosome'),
            value=chromosomes
        )
        workflow.transform(
            name='generate_coverage_bed',
            func='wgs.workflows.postprocessing.tasks.generate_coverage_bed',
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
            func='wgs.workflows.postprocessing.tasks.samtools_coverage',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            axes=('chromosome',),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed', 'chromosome'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),

            ),
            kwargs={'docker_image': config.containers('samtools')}
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


def create_postprocessing_workflow(
        normal_bam,
        tumour_bam,
        titan,
        remixt,
        breakpoints_consensus,
        roh,
        germline_calls,
        somatic_calls,
        circos_plot_remixt,
        circos_plot_titan,
        genome_wide_plot,
        refdir,
        sample_id,
        single_node=False
):

    refdir_paths = config.refdir_data(refdir)['paths']
    refdir_params = config.refdir_data(refdir)['params']

    ideogram = refdir_paths["ideogram"]

    titan_calls = titan[sample_id]
    remixt_calls = remixt[sample_id]
    sv_calls = breakpoints_consensus[sample_id]
    roh_calls = roh[sample_id]
    germline_vcf = germline_calls[sample_id]
    somatic_calls = somatic_calls[sample_id]
    chromosomes = refdir_params['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name='coverage_normal_data',
        func=get_coverage_data,
        args=(
            mgd.InputFile(normal_bam),
            mgd.TempOutputFile('normal_coverage'),
            refdir,
        ),
        kwargs={'single_node': single_node}
    )

    workflow.subworkflow(
        name='coverage_tumour_data',
        func=get_coverage_data,
        args=(
            mgd.InputFile(tumour_bam),
            mgd.TempOutputFile('tumour_coverage'),
            refdir,
        ),
        kwargs={'single_node': single_node}
    )


    workflow.transform(
        name='parse_roh',
        ctx=helpers.get_default_ctx(
            memory=5
        ),
        func="wgs.workflows.postprocessing.tasks.parse_roh",
        args=(
            mgd.InputFile(roh_calls),
            mgd.TempOutputFile("ROH_parsed"),
        ),
    )

    if remixt_calls:

        workflow.transform(
            name='generate_genome_wide_plot',
            ctx=helpers.get_default_ctx(
                memory=10,
            ),
            func="wgs.workflows.postprocessing.tasks.genome_wide",
            args=(
                mgd.InputFile(titan_calls),
                mgd.TempInputFile("ROH_parsed"),
                mgd.InputFile(germline_vcf),
                mgd.InputFile(somatic_calls),
                mgd.TempInputFile('tumour_coverage'),
                mgd.TempInputFile('normal_coverage'),
                mgd.InputFile(sv_calls),
                mgd.InputFile(ideogram),
                chromosomes,
                mgd.OutputFile(genome_wide_plot),

            ),
            kwargs={"remixt": mgd.InputFile(remixt_calls),
                    "remixt_label": sample_id}
        )
        workflow.transform(
            name='generate_circos_plot',
            ctx=helpers.get_default_ctx(
                memory=10
            ),
            func="wgs.workflows.postprocessing.tasks.circos",
            args=(
                mgd.InputFile(titan_calls),
                sample_id,
                mgd.InputFile(sv_calls),
                mgd.TempOutputFile(circos_plot_remixt),
                mgd.TempOutputFile(circos_plot_titan),
                mgd.TempSpace('circos'),
            ),

            kwargs={'docker_image': config.containers('circos'),
                    'remixt_calls': mgd.InputFile(remixt_calls)},
        )
    else:

        workflow.transform(
            name='generate_genome_wide_plot',
            ctx=helpers.get_default_ctx(
                memory=10,
            ),
            func="wgs.workflows.postprocessing.tasks.genome_wide",
            args=(
                mgd.InputFile(titan_calls),
                mgd.TempInputFile("ROH_parsed"),
                mgd.InputFile(germline_vcf),
                mgd.InputFile(somatic_calls),
                mgd.InputFile('tumour_coverage'),
                mgd.InputFile('normal_coverage'),
                mgd.InputFile(sv_calls),
                mgd.InputFile(ideogram),
                chromosomes,
                mgd.OutputFile(genome_wide_plot),
            ),
        )

        workflow.transform(
            name='generate_circos_plot',
            ctx=helpers.get_default_ctx(
                memory=10
            ),
            func="wgs.workflows.postprocessing.tasks.circos",
            args=(
                mgd.InputFile(titan_calls),
                sample_id,
                mgd.InputFile(sv_calls),
                mgd.OutputFile(circos_plot_remixt),
                mgd.OutputFile(circos_plot_titan),
                mgd.TempSpace('circos'),
            ),

            kwargs={'docker_image': config.containers('circos')}
        )

    return workflow

