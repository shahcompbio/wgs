import pypeliner
import pypeliner.managed as mgd
from pypeliner.workflow import Workflow
from wgs.config import config
from wgs.utils import helpers


def create_strelka_workflow(
        normal_bam_file,
        tumour_bam_file,
        snv_vcf_file,
        snv_maf_file,
        indel_vcf_file,
        indel_maf_file,
        reference,
        reference_vep,
        chromosomes,
        normal_id,
        tumour_id,
        single_node=False,
        is_exome=False
):
    params = config.default_params('variant_calling')

    workflow = Workflow(
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00'),
    )

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.mutationseq.tasks.generate_intervals',
        ret=mgd.OutputChunks('regions'),
        args=(
            reference,
            chromosomes
        ),
        kwargs={'size': params['split_size']}
    )

    workflow.transform(
        name='count_fasta_bases',
        func="wgs.workflows.strelka.tasks.count_fasta_bases",
        args=(
            reference,
            pypeliner.managed.TempOutputFile('ref_base_counts.tsv'),
        ),
        kwargs={'docker_image': config.containers('strelka')}
    )

    workflow.transform(
        name="get_chrom_sizes",
        func="wgs.workflows.strelka.tasks.get_known_chromosome_sizes",
        ret=pypeliner.managed.TempOutputObj('known_sizes'),
        args=(
            pypeliner.managed.TempInputFile('ref_base_counts.tsv'),
            chromosomes
        )
    )

    if single_node:
        workflow.transform(
            name='strelka_one_node',
            func="wgs.workflows.strelka.tasks.strelka_one_node",
            args=(
                pypeliner.managed.InputFile(normal_bam_file, extensions=['.bai']),
                pypeliner.managed.InputFile(tumour_bam_file, extensions=['.bai']),
                reference,
                mgd.TempOutputFile('indels.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempOutputFile('snvs.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempSpace('call_genome_segment_tmp'),
                mgd.InputChunks('regions'),
                mgd.TempInputObj('known_sizes'),
            ),
            kwargs={
                'is_exome': is_exome,
                'strelka_docker_image': config.containers('strelka'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name='get_chromosome_depths',
            axes=('regions',),
            func="wgs.workflows.strelka.tasks.get_chromosome_depth",
            args=(
                mgd.InputInstance('regions'),
                pypeliner.managed.InputFile(normal_bam_file, extensions=['.bai']),
                reference,
                mgd.TempOutputFile('chrom_depth.txt', 'regions'),
            ),
            kwargs={'docker_image': config.containers('strelka')},
        )

        workflow.transform(
            name='merge_chromosome_depths',
            func="wgs.workflows.strelka.tasks.merge_chromosome_depths",
            args=(
                mgd.TempInputFile('chrom_depth.txt', 'regions', axes_origin=[]),
                mgd.TempOutputFile('merged_chrom_depth.txt')
            )
        )

        workflow.transform(
            name='call_genome_segment',
            axes=('regions',),
            func="wgs.workflows.strelka.tasks.call_genome_segment",
            args=(
                mgd.TempInputFile('merged_chrom_depth.txt'),
                pypeliner.managed.InputFile(normal_bam_file, extensions=['.bai']),
                pypeliner.managed.InputFile(tumour_bam_file, extensions=['.bai']),
                reference,
                mgd.TempOutputFile('indels.vcf', 'regions'),
                mgd.TempOutputFile('snvs.vcf', 'regions'),
                mgd.TempSpace('call_genome_segment_tmp', 'regions'),
                mgd.InputInstance('regions'),
                mgd.TempInputObj('known_sizes'),
            ),
            kwargs={
                'is_exome': False,
                'docker_image': config.containers('strelka')
            }
        )

        workflow.transform(
            name='merge_indels',
            func='wgs.workflows.strelka.tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('indels.vcf', 'regions'),
                mgd.TempOutputFile('indels.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempSpace("indels_merge")
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

        workflow.transform(
            name='merge_snvs',
            func='wgs.workflows.strelka.tasks.concatenate_vcf',
            args=(
                mgd.TempInputFile('snvs.vcf', 'regions'),
                mgd.TempOutputFile('snvs.vcf.gz', extensions=['.tbi', '.csi']),
                mgd.TempSpace("snvs_merge")
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

    workflow.transform(
        name='filter_vcf_indel',
        func='wgs.workflows.strelka.tasks.filter_vcf',
        args=(
            mgd.TempInputFile('indels.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.transform(
        name='filter_vcf_snv',
        func='wgs.workflows.strelka.tasks.filter_vcf',
        args=(
            mgd.TempInputFile('snvs.vcf.gz', extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.subworkflow(
        name="strelka_snv_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(snv_vcf_file, extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_maf_file),
            reference_vep,
        ),
        kwargs={'tumour_id': tumour_id, 'normal_id': normal_id}
    )

    workflow.subworkflow(
        name="strelka_indel_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(indel_vcf_file, extensions=['.tbi', '.csi']),
            mgd.OutputFile(indel_maf_file),
            reference_vep,
        ),
        kwargs={'tumour_id': tumour_id, 'normal_id': normal_id}
    )

    return workflow
