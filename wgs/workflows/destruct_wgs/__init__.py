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

    workflow.transform(
        name='filter_annotate_breakpoints',
        ctx=helpers.get_default_ctx(memory=4),
        func='wgs.workflows.destruct_wgs.filter_annotate.filter_annotate_breakpoints',
        args=(
            mgd.TempInputFile("raw_breakpoints"),
            mgd.TempInputFile("raw_library"),
            [sample_id + 'N'],
            mgd.TempOutputFile("filter_annotate_breakpoints_output"),
            mgd.TempOutputFile("library"),
        )
    )

    workflow.transform(
        name='mappability_annotate_breakpoints',
        ctx=helpers.get_default_ctx(memory=4),
        func='wgs.workflows.destruct_wgs.flag_mappability.main',
        args=(
            mgd.TempInputFile("filter_annotate_breakpoints_output"),
            mgd.TempOutputFile("breakpoints"),
            sv_config["mappability_ref"],
        )
    )

    workflow.transform(
        name='prep_raw_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("raw_breakpoints"),
            mgd.TempOutputFile("raw_breakpoints.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_raw_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_breakpoints.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(raw_breakpoints, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_raw_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("raw_library"),
            mgd.TempOutputFile("raw_library.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_raw_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_library.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(raw_library, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_reads',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("reads"),
            mgd.TempOutputFile("reads.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_reads',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("reads.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(reads, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("breakpoints"),
            mgd.TempOutputFile("breakpoints.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_breakpoints',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("breakpoints.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(breakpoints, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='prep_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.prep_csv_files",
        args=(
            mgd.TempInputFile("library"),
            mgd.TempOutputFile("library.csv.gz", extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='finalize_library',
        ctx={'mem': 8, 'ncpus': 1},
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("library.csv.gz", extensions=['.yaml']),
            mgd.OutputFile(library, extensions=['.yaml']),
        )
    )

    return workflow
