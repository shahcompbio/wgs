'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_museq_workflow(
        snv_vcf,
        snv_maf,
        museqportrait_pdf,
        reference,
        reference_vep,
        chromosomes,
        sample_id,
        thousand_genomes=None,
        dbsnp=None,
        germline_refdata=None,
        tumour_bam=None,
        normal_bam=None,
        single_node=None
):
    name = 'run_museq'
    if tumour_bam:
        tumour_bam = mgd.InputFile(tumour_bam, extensions=['.bai'])
        name += '_tumour'
    if normal_bam:
        normal_bam = mgd.InputFile(normal_bam, extensions=['.bai'])
        name += '_normal'
    single = False if name == 'run_museq_tumour_normal' else True

    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config.containers('wgs')})

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.mutationseq.tasks.generate_intervals',
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
            name=name,
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus='8',
                disk=600
            ),
            func='wgs.utils.museq_utils.run_museq_one_job',
            args=(
                mgd.TempSpace("run_museq_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                params['museq_params'],
            ),
            kwargs={
                'tumour_bam': tumour_bam,
                'normal_bam': normal_bam,
                'museq_docker_image': config.containers('mutationseq'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name=name,
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.utils.museq_utils.run_museq',
            args=(
                mgd.TempOutputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('museq.log', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                params['museq_params'],
                mgd.TempSpace('museq_temp', 'interval')
            ),
            kwargs={
                'tumour_bam': tumour_bam,
                'normal_bam': normal_bam,
                'docker_image': config.containers('mutationseq'),
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
                mgd.TempInputFile('museq.vcf', 'interval'),
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
            mgd.OutputFile(snv_vcf, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': config.containers('vcftools')}
    )

    workflow.transform(
        name='run_museqportrait',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='8:00',
        ),
        func='wgs.workflows.mutationseq.tasks.run_museqportrait',
        args=(
            mgd.InputFile(snv_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(museqportrait_pdf),
            mgd.TempOutputFile('museqportrait.txt'),
            mgd.TempOutputFile('museqportrait.log'),
            single,
        ),
        kwargs={'docker_image': config.containers('mutationseq'),
                'thousand_genomes': thousand_genomes,
                'dbsnp': dbsnp,
                'germline_refdata': germline_refdata,
                'germline_plot_threshold': params['germline_portrait_threshold']
                }
    )

    if single:
        vcf2maf_kwargs = {'normal_id': sample_id}
    else:
        vcf2maf_kwargs = {'tumour_id': sample_id}
    workflow.subworkflow(
        name="mutationseq_single_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(snv_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_maf),
            reference_vep
        ),
        kwargs=vcf2maf_kwargs
    )

    return workflow
