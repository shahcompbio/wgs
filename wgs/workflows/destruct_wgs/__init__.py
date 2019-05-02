'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd

import tasks
import filter_annotate


def create_destruct_workflow(
        tumour_bam, normal_bam, raw_breakpoints, raw_library,
        breakpoints, library, reads,
        sample_id, global_config, sv_config,
        single_node=False
):

    destruct_config = {}

    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='destruct_local',
            func=tasks.run_destruct_local,
            args=(
                mgd.TempSpace("destruct_local_temp"),
                mgd.InputFile(tumour_bam),
                mgd.InputFile(normal_bam),
                sample_id,
                mgd.OutputFile(raw_breakpoints),
                mgd.OutputFile(raw_library),
                mgd.OutputFile(reads),
                destruct_config,
                sv_config['refdata_destruct'],
            ),
            kwargs={'ncpus': None}
        )
    else:
        workflow.subworkflow(
            name='destruct_parallel',
            func='destruct.workflow.create_destruct_workflow',
            args=(
                {sample_id: mgd.InputFile(tumour_bam),
                 sample_id + 'N': mgd.InputFile(normal_bam)},
                mgd.OutputFile(raw_breakpoints),
                mgd.OutputFile(raw_library),
                mgd.OutputFile(reads),
                destruct_config,
                sv_config['refdata_destruct']
            )
        )

    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx={
            'mem': 4,
            'ncpus': 1,
        },
        func=filter_annotate.filter_annotate_breakpoints,
        args=(
            mgd.InputFile(raw_breakpoints),
            mgd.InputFile(raw_library),
            [sample_id + 'N'],  # control_ids
            mgd.OutputFile(breakpoints),
            mgd.OutputFile(library),
        )
    )

    return workflow
