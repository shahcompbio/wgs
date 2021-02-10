import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


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
                memory=5,
                walltime='24:00'),
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
                memory=5,
                walltime='24:00'),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                mapping_qual,

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
            func='wgs.workflows.sample_qc.tasks.generate_coverage_bed',
            ctx=helpers.get_default_ctx(
                memory=5,
                walltime='24:00'),
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
                memory=5,
                walltime='24:00'),
            axes=('chromosome',),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempInputFile('coverage_bed.bed', 'chromosome'),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                mapping_qual,
            ),
            kwargs={'docker_image': config.containers('samtools')}
        )

        workflow.transform(
            name='merge_data',
            func='wgs.utils.csvutils.concatenate_csv',
            ctx=helpers.get_default_ctx(
                memory=5,
                walltime='6:00'),
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
        name='parse_roh',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='6:00'),
        func="wgs.workflows.sample_qc.tasks.parse_roh",
        args=(
            mgd.InputFile(roh),
            mgd.TempOutputFile("ROH_parsed"),
        ),
    )

    workflow.transform(
        name='generate_genome_wide_plot',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='24:00'),
        func="wgs.workflows.sample_qc.tasks.genome_wide",
        args=(
            sample_id,
            mgd.InputFile(titan),
            mgd.TempInputFile("ROH_parsed"),
            mgd.InputFile(germline_calls),
            mgd.InputFile(somatic_calls),
            mgd.InputFile(remixt),
            mgd.InputFile(tumour_coverage),
            mgd.InputFile(normal_coverage),
            mgd.InputFile(breakpoints_consensus),
            chromosomes,
            mgd.OutputFile(genome_wide_plot),
        ),
    )

    return workflow
