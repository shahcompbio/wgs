'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers

def create_museq_workflow(
        snv_vcf,
        museqportrait_pdf,
        global_config,
        varcall_config,
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

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': varcall_config['docker']['wgs']})

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.mutationseq.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='1:00',
        ),
        ret=mgd.OutputChunks('interval'),
        args=(
            varcall_config['reference'],
            varcall_config['chromosomes']
        ),
        kwargs={'size': varcall_config['split_size']}
    )

    if single_node:
        workflow.transform(
            name=name,
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='48:00',
                ncpus=global_config['threads'],
                disk=600
        ),
            func='wgs.utils.museq_utils.run_museq_one_job',
            args=(
                mgd.TempSpace("run_museq_temp"),
                mgd.TempOutputFile('merged.vcf'),
                varcall_config['reference'],
                mgd.InputChunks('interval'),
                varcall_config['museq_params'],
            ),
            kwargs={
                'tumour_bam': tumour_bam,
                'normal_bam': normal_bam,
                'museq_docker_image': varcall_config['docker']['mutationseq'],
                'vcftools_docker_image': varcall_config['docker']['vcftools']
            }
        )
    else:
        workflow.transform(
            name=name,
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.utils.museq_utils.run_museq',
            args=(
                mgd.TempOutputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('museq.log', 'interval'),
                varcall_config['reference'],
                mgd.InputInstance('interval'),
                varcall_config['museq_params'],
            ),
            kwargs={
                'tumour_bam': tumour_bam,
                'normal_bam': normal_bam,
                'docker_image': varcall_config['docker']['mutationseq']
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='8:00',
            ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            ),
            kwargs={'docker_image': varcall_config['docker']['vcftools']}
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
        kwargs={'docker_image': varcall_config['docker']['vcftools']}
    )

    workflow.transform(
        name='run_museqportrait',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='8:00',
        ),
        func='wgs.workflows.mutationseq.tasks.run_museqportrait',
        args=(
            mgd.InputFile(snv_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(museqportrait_pdf),
            mgd.TempOutputFile('museqportrait.txt'),
            mgd.TempOutputFile('museqportrait.log'),
            single,
            varcall_config['plot_params'],
        ),
        kwargs={'docker_image': varcall_config['docker']['mutationseq']}
    )

    return workflow
