'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import vcf_tasks

import tasks


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

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='generate_intervals',
        func=tasks.generate_intervals,
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '01:00'},
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
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '24:00'},
            func=tasks.run_museq_one_job,
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
            }
        )
    else:
        workflow.transform(
            name=name,
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '24:00'},
            axes=('interval',),
            func=tasks.run_museq,
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
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '08:00'},
            func=tasks.merge_vcfs,
            args=(
                mgd.TempInputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            )
        )

    workflow.transform(
        name='finalise_snvs',
        ctx={'ncpus': 1, 'walltime': '01:00'},
        func=vcf_tasks.finalise_vcf,
        args=(
            mgd.TempInputFile('merged.vcf'),
            mgd.OutputFile(snv_vcf, extensions=['.tbi', '.csi']),
        ),
    )

    workflow.transform(
        name='run_museqportrait',
        ctx={'num_retry': 3, 'mem_retry_increment': 2,
             'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '08:00'},
        func=tasks.run_museqportrait,
        args=(
            mgd.InputFile(snv_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(museqportrait_pdf),
            mgd.TempOutputFile('museqportrait.txt'),
            mgd.TempOutputFile('museqportrait.log'),
            single,
            varcall_config['plot_params'],
        ),
    )

    return workflow
