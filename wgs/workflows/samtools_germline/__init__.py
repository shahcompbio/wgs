'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers

def create_samtools_germline_workflow(
        germline_vcf,
        global_config,
        varcall_config,
        bam_file,
        single_node=None
):

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': varcall_config['docker']['wgs']})

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.samtools_germline.tasks.generate_intervals',
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
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='48:00',
                ncpus=global_config['threads'],
                disk=600
        ),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline_one_job',
            args=(
                mgd.TempSpace("run_samtools_temp"),
                mgd.TempOutputFile('merged.vcf'),
                varcall_config['reference'],
                mgd.InputChunks('interval'),
                mgd.InputFile(bam_file)
            ),
            kwargs={
                'docker_image': varcall_config['docker']['samtools']
            }
        )
    else:
        workflow.transform(
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline',
            args=(
                mgd.TempOutputFile('germline.vcf', 'interval'),
                varcall_config['reference'],
                mgd.InputInstance('interval'),
                mgd.InputFile(bam_file)
            ),
            kwargs={
                'docker_image': varcall_config['docker']['samtools']
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
                mgd.TempInputFile('germline.vcf', 'interval'),
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
            mgd.OutputFile(germline_vcf, extensions=['.tbi', '.csi']),
        ),
        kwargs={'docker_image': varcall_config['docker']['vcftools']}
    )

    return workflow
