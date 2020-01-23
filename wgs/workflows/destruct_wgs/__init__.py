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

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': sv_config['docker']['wgs']})

    if single_node:
        workflow.transform(
            name='destruct_local',
            ctx=helpers.get_default_ctx(
                walltime='120:00',
                disk=800
            ),
            func='wgs.workflows.destruct_wgs.tasks.run_destruct_local',
            args=(
                mgd.TempSpace("destruct_local_temp"),
                mgd.InputFile(tumour_bam),
                mgd.InputFile(normal_bam),
                sample_id,
                mgd.TempOutputFile("raw_breakpoints"),
                mgd.TempOutputFile("raw_library"),
                mgd.TempOutputFile("reads"),
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
                walltime="48:00",
            ),
            # refers to seperate destruct package
            func='destruct.workflow.create_destruct_workflow',
            args=(
                {sample_id: mgd.InputFile(tumour_bam),
                 sample_id + 'N': mgd.InputFile(normal_bam)},
                mgd.TempOutputFile("raw_breakpoints"),
                mgd.TempOutputFile("raw_library"),
                mgd.TempOutputFile("reads"),
                destruct_config,
                sv_config['refdata_destruct']
            )
        )

    workflow.commandline(
        name='filter_annotate_breakpoints',
        ctx=helpers.get_default_ctx(
            docker_image=sv_config['docker']['destruct'],
            memory=4,
            walltime='8:00'
        ),
        args=(
            'filter_annotate_breakpoints.py',
            '--breakpoints',
            mgd.TempInputFile("raw_breakpoints"),
            '--library',
            mgd.TempInputFile("raw_library"),
            '--control_ids',
            sample_id + 'N',
            '--out_breakpoints',
            mgd.TempOutputFile("filter_annotate_breakpoints_output"),
            '--out_library',
            mgd.TempOutputFile("library"),
        )
    )

    workflow.transform(
        name='mappability_annotate_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func='wgs.workflows.destruct_wgs.flag_mappability.main',
        args=(
            mgd.TempInputFile("filter_annotate_breakpoints_output"),
            mgd.TempOutputFile("breakpoints"),
            sv_config["mappability_ref"],
        )
    )

    workflow.transform(
        name='prep_raw_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("raw_breakpoints"),
            mgd.TempOutputFile("raw_breakpoints.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_raw_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_breakpoints.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(raw_breakpoints, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_raw_library',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("raw_library"),
            mgd.TempOutputFile("raw_library.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_raw_library',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_library.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(raw_library, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_reads',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("reads"),
            mgd.TempOutputFile("reads.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_reads',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(reads, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("breakpoints"),
            mgd.TempOutputFile("breakpoints.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("breakpoints.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(breakpoints, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_library',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("library"),
            mgd.TempOutputFile("library.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_library',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("library.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(library, extensions=['.yaml']),
        )
    )

    return workflow
