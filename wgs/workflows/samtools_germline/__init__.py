'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_samtools_germline_workflow(
        germline_vcf,
        germline_maf,
        germline_roh,
        bam_file,
        reference,
        reference_vep,
        chromosomes,
        sample_id,
        single_node=None
):
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config.containers('wgs')})

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.samtools_germline.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='1:00',
        ),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            chromosomes
        ),
        kwargs={'size': params['split_size']}
    )

    if single_node:
        workflow.transform(
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline_one_job',
            args=(
                mgd.TempSpace("run_samtools_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(bam_file)
            ),
            kwargs={
                'samtools_docker_image': config.containers('samtools'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline',
            args=(
                mgd.TempOutputFile('germline.vcf.gz', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(bam_file),
                mgd.TempSpace('tempdir_samtools', 'interval')
            ),
            kwargs={
                'docker_image': config.containers('samtools')
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='8:00',
            ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('germline.vcf.gz', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

    workflow.transform(
        name='finalise_snvs',
        ctx=helpers.get_default_ctx(
            walltime='8:00',
        ),
        func='wgs.utils.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('merged.vcf'),
            mgd.OutputFile(germline_vcf, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.transform(
        name='roh_calling',
        ctx=helpers.get_default_ctx(
            walltime='8:00',
        ),
        func='wgs.workflows.samtools_germline.tasks.roh_calling',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_roh)
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.subworkflow(
        name="samtools_germline_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_maf, extensions=['.tbi', '.csi']),
            reference_vep
        ),
        kwargs={'normal_id': sample_id}
    )


    return workflow
