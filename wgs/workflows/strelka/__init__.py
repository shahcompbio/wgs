import pypeliner.managed as mgd
from pypeliner.workflow import Workflow
from wgs.utils import helpers

from wgs.config import config

def create_strelka_workflow(
        normal_bam,
        tumour_bam,
        indel_vcf_file,
        snv_vcf_file,
        reference,
        chromosomes,
        use_depth_thresholds=True,
        single_node=False,
):

    params = config.default_params('variant_calling')

    workflow = Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.strelka.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00'),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            chromosomes
        ),
        kwargs={'size': params['split_size']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.transform(
        name='count_fasta_bases',
        ctx={
            'mem': '5',
            'walltime': '4:00'
        },
        func='wgs.workflows.strelka.tasks.count_fasta_bases',
        args=(
            reference,
            mgd.TempOutputFile('ref_base_counts.tsv'),
        ),
        kwargs={'docker_image': config.containers('strelka')}
    )

    workflow.transform(
        name='get_chrom_sizes',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00'),
        func='wgs.workflows.strelka.tasks.get_known_chromosome_sizes',
        ret=mgd.TempOutputObj('known_sizes'),
        args=(
            mgd.TempInputFile('ref_base_counts.tsv'),
        ),
    )

    if single_node:
        workflow.transform(
            name='call_somatic_variants',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='96:00',
                ncpus=8),
            func='wgs.workflows.strelka.tasks.call_somatic_variants_one_job',
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
                mgd.TempInputObj('known_sizes'),
                reference,
                mgd.InputChunks('interval'),
                chromosomes,
                mgd.TempSpace("strelka_single_node_run"),
            ),
            kwargs={
                'strelka_docker_image': config.containers('strelka'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name='call_somatic_variants',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='16:00', ),
            axes=('interval',),
            func='wgs.workflows.strelka.tasks.call_somatic_variants',
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempInputObj('known_sizes'),
                reference,
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf.window', 'interval'),
                mgd.TempOutputFile('somatic.snvs.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('strelka.stats', 'interval'),
                mgd.InputInstance('interval'),
            ),
            kwargs={'docker_image': config.containers('strelka')}
        )

        workflow.transform(
            name='add_indel_filters',
            axes=('chrom',),
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='8:00', ),
            func='wgs.workflows.strelka.tasks.filter_indel_file_list',
            args=(
                mgd.TempInputFile('somatic.indels.unfiltered.vcf', 'interval', axes_origin=[]),
                mgd.TempInputFile('strelka.stats', 'interval', axes_origin=[]),
                mgd.TempInputFile('somatic.indels.unfiltered.vcf.window', 'interval', axes_origin=[]),
                mgd.TempOutputFile('somatic.indels.filtered.vcf', 'chrom'),
                mgd.InputInstance('chrom'),
                mgd.TempInputObj('known_sizes'),
                mgd.InputChunks('interval'),
            ),
            kwargs={'use_depth_filter': use_depth_thresholds}
        )

        workflow.transform(
            name='add_snv_filters',
            axes=('chrom',),
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='8:00', ),
            func='wgs.workflows.strelka.tasks.filter_snv_file_list',
            args=(
                mgd.TempInputFile('somatic.snvs.unfiltered.vcf', 'interval', axes_origin=[]),
                mgd.TempInputFile('strelka.stats', 'interval', axes_origin=[]),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf', 'chrom'),
                mgd.InputInstance('chrom'),
                mgd.TempInputObj('known_sizes'),
                mgd.InputChunks('interval'),
            ),
            kwargs={'use_depth_filter': use_depth_thresholds}
        )

        workflow.transform(
            name='merge_indels',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='8:00', ),
            func='wgs.workflows.strelka.vcf_tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('somatic.indels.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
                mgd.TempSpace('merge_indels_tempdir')
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

        workflow.transform(
            name='merge_snvs',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='8:00', ),
            func='wgs.workflows.strelka.vcf_tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('somatic.snvs.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
                mgd.TempSpace('merge_snvs_tempdir')
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

    workflow.transform(
        name='filter_indels',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        func='wgs.workflows.strelka.vcf_tasks.filter_vcf',
        args=(
            mgd.TempInputFile('somatic.indels.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.indels.passed.vcf'),
        ),

    )

    workflow.transform(
        name='filter_snvs',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        func='wgs.workflows.strelka.vcf_tasks.filter_vcf',
        args=(
            mgd.TempInputFile('somatic.snvs.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.snvs.passed.vcf'),
        ),
    )

    workflow.transform(
        name='finalise_indels',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        func='wgs.workflows.strelka.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('somatic.indels.passed.vcf'),
            mgd.OutputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.transform(
        name='finalise_snvs',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        func='wgs.workflows.strelka.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('somatic.snvs.passed.vcf'),
            mgd.OutputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    return workflow
