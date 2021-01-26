'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_destruct_wgs_workflow(
        tumour_bam, normal_bam, raw_breakpoints, raw_library,
        breakpoints, library, reads,
        sample_id, reference, destruct_refdata, gtf, mappability,
        single_node=False
):
    destruct_config = {
        'genome_fasta': reference,
        'genome_fai': reference + '.fai',
        'gtf_filename': gtf
    }

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': config.containers('wgs')})

    workflow.transform(
        name="get_destruct_config",
        func="destruct.defaultconfig.get_config",
        ctx=helpers.get_default_ctx(
            docker_image=config.containers('destruct'),
            walltime="48:00",
        ),
        ret=mgd.TempOutputObj("destruct_config"),
        args=(
            destruct_refdata,
            destruct_config
        )
    )

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
                mgd.TempOutputFile("raw_reads"),
                mgd.TempInputObj("destruct_config"),
                destruct_refdata,
            ),
            kwargs={'ncpus': 16, 'docker_image': config.containers('destruct')}
        )
    else:
        workflow.subworkflow(
            name='destruct_parallel',
            ctx=helpers.get_default_ctx(
                docker_image=config.containers('destruct'),
                walltime="48:00",
            ),
            # refers to seperate destruct package
            func='destruct.workflow.create_destruct_workflow',
            args=(
                {sample_id: mgd.InputFile(tumour_bam),
                 sample_id + 'N': mgd.InputFile(normal_bam)},
                mgd.TempOutputFile("raw_breakpoints"),
                mgd.TempOutputFile("raw_library"),
                mgd.TempOutputFile("raw_reads"),
                mgd.TempInputObj("destruct_config"),
                destruct_refdata
            )
        )

    workflow.commandline(
        name='filter_annotate_breakpoints',
        ctx=helpers.get_default_ctx(
            docker_image=config.containers('destruct'),
            memory=8,
            walltime='8:00'
        ),
        args=(
            'destruct',
            'extract_somatic',
            mgd.TempInputFile("raw_breakpoints"),
            mgd.TempInputFile("raw_library"),
            mgd.TempOutputFile("filter_annotate_breakpoints_output"),
            mgd.TempOutputFile("library"),
            '--control_ids',
            sample_id + 'N',
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
            mappability,
        )
    )

    workflow.transform(
        name='finalize_raw_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_breakpoints"),
            mgd.OutputFile(raw_breakpoints, extensions=['.yaml']),
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
            mgd.TempInputFile("raw_library"),
            mgd.OutputFile(raw_library, extensions=['.yaml']),
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
            mgd.TempInputFile("breakpoints"),
            mgd.OutputFile(breakpoints, extensions=['.yaml']),
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
            mgd.TempInputFile("library"),
            mgd.OutputFile(library, extensions=['.yaml']),
        )
    )

    workflow.transform(
        name='reheader_reads',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.workflows.destruct_wgs.tasks.reheader_reads",
        args=(
            mgd.TempInputFile("raw_reads"),
            mgd.TempOutputFile("raw_reads_reheader"),
        )
    )

    workflow.transform(
        name='finalize_reads',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00'
        ),
        func="wgs.utils.csvutils.finalize_csv",
        args=(
            mgd.TempInputFile("raw_reads_reheader"),
            mgd.OutputFile(reads, extensions=['.yaml']),
        )
    )

    return workflow
