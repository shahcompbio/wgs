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
        sample_id):
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
            mgd.TempOutputFile('lumpy_filt.csv'),
            svcalling_config["parse_lumpy"],
            sample_id
        ),
        kwargs={'docker_image': svcalling_config['docker']['vizutils']}
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
            mgd.TempInputFile('lumpy.csv'),
            mgd.TempOutputFile('destruct.csv'),
            mgd.OutputFile(output),
            svcalling_config["parse_destruct"],
            sample_id
        ),
        kwargs={'docker_image': svcalling_config['docker']['vizutils']}
    )

    return workflow
