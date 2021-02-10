'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_freebayes_germline_workflow(
        germline_vcf,
        germline_maf,
        bam_file,
        reference,
        reference_vep,
        chromosomes,
        normal_id,
        single_node=None
):
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config.containers('wgs')})

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.freebayes.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='6:00',
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
            name='freebayes_one_node',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            func='wgs.workflows.freebayes.tasks.run_freebayes_one_job',
            args=(
                mgd.TempSpace("run_freebayes_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(bam_file)
            ),
            kwargs={
                'freebayes_docker_image': config.containers('freebayes'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name='freebayes',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.freebayes.tasks.run_freebayes_germline',
            args=(
                mgd.TempOutputFile('freebayes_germline.vcf', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(bam_file),
                mgd.TempSpace('tempdir_freebayes', 'interval')
            ),
            kwargs={
                'docker_image': config.containers('freebayes')
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='6:00',
            ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('freebayes_germline.vcf', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

    workflow.transform(
        name='finalise_snvs',
        ctx=helpers.get_default_ctx(
            walltime='6:00',
        ),
        func='wgs.utils.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('merged.vcf'),
            mgd.OutputFile(germline_vcf, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.subworkflow(
        name="freebayes_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_maf, extensions=['.tbi', '.csi']),
            reference_vep,
        ),
        kwargs={'normal_id': normal_id}
    )

    return workflow
