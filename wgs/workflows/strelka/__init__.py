import pypeliner.managed as mgd
from pypeliner.workflow import Workflow
from wgs.utils import helpers

default_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']


def create_strelka_workflow(
        normal_bam,
        tumour_bam,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        global_config,
        varcall_config,
        chromosomes=default_chromosomes,
        use_depth_thresholds=True,
        single_node=False,
):
    docker_config = varcall_config['docker']

    workflow = Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.strelka.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='2:00'),
        ret=mgd.OutputChunks('interval'),
        args=(
            varcall_config['reference'],
            varcall_config['chromosomes']
        ),
        kwargs={'size': varcall_config['split_size']}
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=varcall_config['chromosomes'],
    )

    workflow.transform(
        name='count_fasta_bases',
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '01:00'},
        func='wgs.workflows.strelka.tasks.count_fasta_bases',
        args=(
            ref_genome_fasta_file,
            mgd.TempOutputFile('ref_base_counts.tsv'),
        ),
        kwargs={'docker_image': docker_config['strelka']}
    )

    workflow.transform(
        name='get_chrom_sizes',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='2:00'),
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
                memory=global_config['memory']['med'],
                walltime='48:00',
                ncpus=global_config['threads']),
            func='wgs.workflows.strelka.tasks.call_somatic_variants_one_job',
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
                mgd.TempInputObj('known_sizes'),
                ref_genome_fasta_file,
                mgd.InputChunks('interval'),
                chromosomes,
                mgd.TempSpace("strelka_single_node_run"),
            ),
            kwargs={
                'strelka_docker_image': docker_config['strelka'],
                'vcftools_docker_image': docker_config['vcftools']
            }
        )
    else:
        workflow.transform(
            name='call_somatic_variants',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['med'],
                walltime='8:00', ),
            axes=('interval',),
            func='wgs.workflows.strelka.tasks.call_somatic_variants',
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempInputObj('known_sizes'),
                ref_genome_fasta_file,
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf.window', 'interval'),
                mgd.TempOutputFile('somatic.snvs.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('strelka.stats', 'interval'),
                mgd.InputInstance('interval'),
            ),
            kwargs={'docker_image': docker_config['strelka']}
        )

        workflow.transform(
            name='add_indel_filters',
            axes=('chrom',),
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['med'],
                walltime='2:00', ),
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
                memory=global_config['memory']['med'],
                walltime='2:00', ),
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
                memory=global_config['memory']['med'],
                walltime='2:00', ),
            func='wgs.workflows.strelka.vcf_tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('somatic.indels.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
                mgd.TempSpace('merge_indels_tempdir')
            ),
        )

        workflow.transform(
            name='merge_snvs',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['med'],
                walltime='2:00', ),
            func='wgs.workflows.strelka.vcf_tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('somatic.snvs.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
                mgd.TempSpace('merge_snvs_tempdir')
            ),
        )

    workflow.transform(
        name='filter_indels',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='2:00', ),
        func='wgs.workflows.strelka.vcf_tasks.filter_vcf',
        args=(
            mgd.TempInputFile('somatic.indels.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.indels.passed.vcf'),
        ),
    )

    workflow.transform(
        name='filter_snvs',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='2:00', ),
        func='wgs.workflows.strelka.vcf_tasks.filter_vcf',
        args=(
            mgd.TempInputFile('somatic.snvs.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.snvs.passed.vcf'),
        ),
    )

    workflow.transform(
        name='finalise_indels',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='2:00', ),
        func='wgs.workflows.strelka.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('somatic.indels.passed.vcf'),
            mgd.OutputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': docker_config['vcftools']}
    )

    workflow.transform(
        name='finalise_snvs',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='2:00', ),
        func='wgs.workflows.strelka.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('somatic.snvs.passed.vcf'),
            mgd.OutputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': docker_config['vcftools']}
    )

    return workflow
