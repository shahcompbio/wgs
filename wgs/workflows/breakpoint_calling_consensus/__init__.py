'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers

from wgs.config import config

def create_consensus_workflow(
        destruct_breakpoints,
        lumpy_vcf,
        output,
        chromosomes
):

    params = config.default_params('breakpoint_calling')
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_lumpy',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00',
        ),
        func='wgs.workflows.breakpoint_calling_consensus.tasks.parse_lumpy_task',
        args=(
            mgd.InputFile(lumpy_vcf),
            mgd.TempOutputFile('lumpy.csv'),
            params["parse_lumpy"],
        ),
        kwargs={'chromosomes': chromosomes}
    )

    workflow.transform(
        name='parse_destruct',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00',
        ),
        func='wgs.workflows.breakpoint_calling_consensus.tasks.parse_destruct_task',
        args=(
            mgd.InputFile(destruct_breakpoints),
            mgd.TempOutputFile('destruct.csv'),
            params["parse_destruct"],
        ),
        kwargs={'chromosomes': chromosomes}
    )

    workflow.transform(
        name='consensus_breakpoint_calling',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00',
        ),
        func='wgs.workflows.breakpoint_calling_consensus.tasks.consensus_calls',
        args=(
            mgd.TempInputFile('destruct.csv'),
            mgd.TempInputFile('lumpy.csv'),
            mgd.OutputFile(output, extensions=['.yaml']),
            params['consensus']
        ),
    )

    return workflow
