'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
import tasks


def create_consensus_workflow(
        destruct_breakpoints,
        lumpy_vcf,
        output,
        config):

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_lumpy',
        ctx={'num_retry': 3, 'mem_retry_increment': 2,
            'mem': config['memory']['high'],
            'pool_id': config['pools']['highmem'],
            'ncpus': 1,'walltime': '08:00'},
        func=tasks.parse_lumpy,
        args=(
            mgd.InputFile(lumpy_vcf),
            mgd.TempOutputFile('lumpy.csv'),
            mgd.TempOutputFile('lumpy_filt.csv'),
            config["parse_lumpy"],
        ),
    )


    workflow.transform(
        name='parse_destruct',
        ctx={'num_retry': 3, 'mem_retry_increment': 2,
            'mem': config['memory']['high'],
            'pool_id': config['pools']['highmem'],
            'ncpus': 1,'walltime': '08:00'},
        func=tasks.parse_destruct,
        args=(
            mgd.InputFile(destruct_breakpoints),
            mgd.TempInputFile('lumpy.csv'),
            mgd.TempOutputFile('destruct.csv'),
            mgd.OutputFile(output),
            config["parse_destruct"],
        ),
    )

    return workflow
