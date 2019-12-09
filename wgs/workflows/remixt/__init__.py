'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

import tasks
from wgs.utils import helpers

def create_remixt_workflow(
        tumour_path,
        normal_path,
        breakpoints,
        sample_id,
        remixt_refdata,
        remixt_results_filename,
        remixt_raw_dir,
        min_num_reads,
        global_config,
        config,
        single_node=False,
):
    ctx = {'docker_image': config['docker']['wgs']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    remixt_config = {}

    if breakpoints is None:
        workflow.setobj(
            obj=mgd.TempOutputObj('emptybreakpoints'),
            value=[],
        )

        workflow.transform(
            name='write_empty_breakpoints',
            func=tasks.write_empty_breakpoints,
            args=(
                mgd.TempInputObj('emptybreakpoints'),
                mgd.TempOutputFile('filtered_breakpoints.csv'),
            ),
        )

    else:
        workflow.transform(
            name='filter_breakpoints',
            func=tasks.filter_destruct_breakpoints,
            ctx=helpers.get_default_ctx(
                memory=4,
                walltime='2:00'),
            args=(
                mgd.InputFile(breakpoints),
                mgd.TempOutputFile('filtered_breakpoints.csv'),
            )
        )

    if single_node:
        workflow.transform(
            name='remixt',
            func=tasks.run_remixt_local,
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='72:00',
                ncpus=global_config['threads']),
            args=(
                mgd.TempSpace("remixt_temp"),
                mgd.TempInputFile('filtered_breakpoints.csv'),
                mgd.InputFile(tumour_path, extensions=['.bai']),
                mgd.InputFile(normal_path, extensions=['.bai']),
                sample_id,
                mgd.OutputFile(remixt_results_filename),
                remixt_raw_dir,
                remixt_config,
                remixt_refdata,
            ),
        )
    else:
        workflow.subworkflow(
            name='remixt',
            func="remixt.workflow.create_remixt_bam_workflow",
            ctx={'docker_image': config['docker']['remixt']},
            args=(
                mgd.TempInputFile('filtered_breakpoints.csv'),
                {sample_id: mgd.InputFile(tumour_path, extensions=['.bai']),
                 sample_id + 'N': mgd.InputFile(normal_path, extensions=['.bai'])},
                {sample_id: mgd.OutputFile(remixt_results_filename)},
                remixt_raw_dir,
                remixt_config,
                remixt_refdata,
            ),
            kwargs={
                'normal_id': sample_id + 'N',
            }
        )

    return workflow
