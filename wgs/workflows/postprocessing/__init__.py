import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
import os

def get_coverage_data(input_bam, output, refdir, single_node=False):
    refdir_paths = config.refdir_data(refdir)['paths']
    chromosomes = config.refdir_data(refdir)['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='samtools_coverage',
            func='wgs.workflows.postprocessing.tasks.samtools_coverage',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                mgd.InputFile(input_bam),
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                chromosomes,
                refdir_paths['reference'],
            ),
            kwargs={'bins_per_chrom': 2000},
        )

    else:
        workflow.setobj(
            obj=mgd.OutputChunks('chromosome'),
            value=chromosomes
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
                mgd.TempOutputFile('per_interval.txt', 'chromosome'),
                mgd.InputInstance('chromosome'),
                refdir_paths['reference'],
            ),
            kwargs={'bins_per_chrom': 2000},
        )
        workflow.transform(
            name='merge_data',
            func='wgs.utils.csvutils.concatenate_csv',
            ctx=helpers.get_default_ctx(
                memory=5
            ),
            args=(
                mgd.TempInputFile('per_interval.txt', 'chromosome', axes_origin=[]),
                mgd.OutputFile(output)
            )
        )

    return workflow


def create_postprocessing_workflow(
        normal_bam,
        tumour_bam,
        variant_dir,
        copynumber_dir,
        breakpoint_dir,
        circos_plot,
        genome_wide_plot,
        refdir,
        sample_id,
        single_node=False
):
    refdir_paths = config.refdir_data(refdir)['paths']

    variant_dir = variant_dir[sample_id]
    breakpoint_dir = breakpoint_dir[sample_id]
    copynumber_dir = copynumber_dir[sample_id]

    cn_calls = os.path.join(copynumber_dir, 'titan', '{}_titan_markers.csv.gz'.format(sample_id))
    sv_calls = os.path.join(breakpoint_dir, '{}_filtered_consensus_calls.csv.gz'.format(sample_id))
    roh_calls = os.path.join(variant_dir, '{}_roh.csv.gz'.format(sample_id))
    germline_vcf = os.path.join(variant_dir, '{}_samtools_germline.vcf.gz'.format(sample_id))
    snv_calls = os.path.join(variant_dir, '{}_museq_paired_annotated.vcf.gz'.format(sample_id))

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
        name='generate_circos_plot',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='24:00',
            disk=400
        ),
        func="wgs.workflows.postprocessing.tasks.circos",
        args=(
            mgd.InputFile(cn_calls),
            mgd.InputFile(sv_calls),
            mgd.TempOutputFile(circos_plot),
            mgd.TempSpace('circos'),
            refdir_paths['annotation_genes'],
        ),
        kwargs={
            'docker_image': config.containers('circos'),
        }
    )

    workflow.transform(
        name='generate_genome_wide_plot',
        ctx=helpers.get_default_ctx(
            memory=10,
        ),
        func="wgs.workflows.postprocessing.tasks.genome_wide_plot",
        args=(
            mgd.InputFile(cn_calls),
            mgd.InputFile(roh_calls),
            mgd.InputFile(germline_vcf),
            mgd.InputFile(snv_calls),
            mgd.TempInputFile('tumour_coverage'),
            mgd.TempInputFile('tumour_coverage'),
            mgd.OutputFile(genome_wide_plot),
        ),
    )

    return workflow
