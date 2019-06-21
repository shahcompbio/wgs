'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_destruct_wgs_workflow(
        tumour_bam, normal_bam, raw_breakpoints, raw_library,
        breakpoints, library, reads,
        sample_id, global_config, sv_config,
        single_node=False
):

    destruct_config = {}

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': sv_config['docker']['destruct']})

    if single_node:
        workflow.transform(
            name='destruct_local',
            ctx=helpers.get_default_ctx(
                walltime='72:00',
                disk=800
            ),
            func='wgs.workflows.destruct_wgs.tasks.run_destruct_local',
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
            kwargs={'ncpus': None, 'docker_image': sv_config['docker']['destruct']}
        )
    else:
        workflow.subworkflow(
            name='destruct_parallel',
            ctx=helpers.get_default_ctx(
                docker_image=sv_config['docker']['destruct'],
            ),
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
        ctx=helpers.get_default_ctx(memory=4),
        func='wgs.workflows.destruct_wgs.filter_annotate.filter_annotate_breakpoints',
        args=(
            mgd.InputFile(raw_breakpoints),
            mgd.InputFile(raw_library),
            [sample_id + 'N'],
            mgd.OutputFile(breakpoints),
            mgd.OutputFile(library),
        )
    )

    return workflow
