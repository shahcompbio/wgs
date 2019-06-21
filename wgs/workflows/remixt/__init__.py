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
        single_node=False,
        docker_containers={}
):

    ctx = {'docker_image': docker_containers.get('remixt')}
    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    if not breakpoints:
        # just setting a random obj so that we dont return empty workflow
        workflow.setobj(
            obj=mgd.OutputChunks('some_random_value'),
            value='another_random_string')
        return workflow

    remixt_config = {}

    workflow.transform(
        name='filter_breakpoints',
        func=tasks.filter_destruct_breakpoints,
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='2:00'),
        args=(
            mgd.InputFile(breakpoints),
            mgd.TempOutputFile('filtered_breakpoints.csv'),
            min_num_reads,
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
                mgd.InputFile(breakpoints),
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
            args=(
                mgd.InputFile(breakpoints),
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
