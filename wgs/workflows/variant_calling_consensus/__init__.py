'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd

import tasks
from wgs.utils import helpers

def create_consensus_workflow(
        museq_germline,
        museq_snv,
        strelka_snv,
        strelka_indel,
        output,
        global_config,
        varcall_config):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_museq_snv',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.parse_museq,
        args=(
            mgd.InputFile(museq_snv),
            mgd.TempOutputFile('museq_snv.csv'),
            mgd.TempOutputFile('museq_snv_filt.csv'),
            varcall_config["parse_museq"],
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    workflow.transform(
        name='parse_museq_germlines',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.parse_museq,
        args=(
            mgd.InputFile(museq_germline),
            mgd.TempOutputFile('museq_germlines.csv'),
            mgd.TempOutputFile('museq_germlines_filt.csv'),
            varcall_config["parse_museq"],
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    workflow.transform(
        name='parse_strelka_snv',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.parse_strelka,
        args=(
            mgd.InputFile(strelka_snv),
            mgd.TempOutputFile('strelka_snv.csv'),
            mgd.TempOutputFile('strelka_snv_filt.csv'),
            varcall_config["parse_strelka"],
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    workflow.transform(
        name='parse_strelka_indel',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.parse_strelka,
        args=(
            mgd.InputFile(strelka_indel),
            mgd.TempOutputFile('strelka_indel.csv'),
            mgd.TempOutputFile('strelka_indel_filt.csv'),
            varcall_config["parse_strelka"],
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    workflow.transform(
        name='merge_snvs',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.merge_overlap,
        args=(
            [mgd.TempInputFile('strelka_snv.csv'),
             mgd.TempInputFile('museq_snv.csv')],
            mgd.TempOutputFile('merged_snv.csv'),
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    workflow.transform(
        name='concatenate',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.concatenate,
        args=(
            [mgd.TempInputFile('merged_snv.csv'),
             mgd.TempInputFile('museq_germlines.csv'),
             mgd.TempInputFile('strelka_indel.csv'),
             ],
            mgd.OutputFile(output),
        ),
        kwargs={'docker_image': varcall_config['docker']['vizutils']}
    )

    return workflow
