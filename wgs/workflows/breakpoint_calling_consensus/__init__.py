'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd

import tasks
from wgs.utils import helpers


def create_consensus_workflow(
        destruct_breakpoints,
        lumpy_vcf,
        output,
        global_config,
        svcalling_config,
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_lumpy',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00',
        ),
        func=tasks.parse_lumpy,
        args=(
            mgd.InputFile(lumpy_vcf),
            mgd.TempOutputFile('lumpy.csv'),
            svcalling_config["parse_lumpy"],
        ),
    )

    workflow.transform(
        name='parse_destruct',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00',
        ),
        func=tasks.parse_destruct,
        args=(
            mgd.InputFile(destruct_breakpoints),
            mgd.TempOutputFile('destruct.csv'),
            svcalling_config["parse_destruct"],
        ),
    )

    workflow.transform(
        name='consensus_breakpoint_calling',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00',
        ),
        func=tasks.consensus,
        args=(
            mgd.TempInputFile('destruct.csv'),
            mgd.TempInputFile('lumpy.csv'),
            mgd.OutputFile(output),
        ),
    )

    return workflow
