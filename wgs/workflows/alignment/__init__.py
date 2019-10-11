import os

import biowrappers.components
import biowrappers.components.io
import biowrappers.components.io.bam.tasks
import biowrappers.components.io.fastq.tasks
import biowrappers.pipelines.realignment.tasks
import pypeliner
import pypeliner.managed as mgd
import tasks
from wgs.utils import helpers

import biowrappers


def align_samples(
        config,
        config_globals,
        fastqs_r1,
        fastqs_r2,
        bam_outputs,
        outdir,
        single_node=False
):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane_id'),
        value=fastqs_r1.keys(),
    )

    workflow.subworkflow(
        name='align_samples',
        func=align_sample,
        axes=('sample_id', 'lane_id'),
        args=(
            config,
            mgd.InputFile('input.r1.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r2),
            mgd.TempOutputFile('aligned_lanes.bam', 'sample_id', 'lane_id'),
            outdir,
            [mgd.InputInstance("sample_id"),
             mgd.InputInstance("lane_id")]
        ),
        kwargs={'single_node': single_node}
    )

    workflow.transform(
        name='merge_tumour_lanes',
        ctx=helpers.get_default_ctx(
            memory=config_globals['memory']['med'],
            walltime='24:00',
            disk=400
        ),
        func="wgs.workflows.alignment.tasks.merge_bams",
        axes=('sample_id',),
        args=(
            mgd.TempInputFile('aligned_lanes.bam', 'sample_id', 'lane_id'),
            mgd.OutputFile('merged_lanes.bam', 'sample_id', fnames=bam_outputs),
        ),
        kwargs={
            'picard_docker_image': config['docker']['picard'],
            'samtools_docker_image': config['docker']['samtools']
        }
    )

    return workflow


def align_sample(
        config, fastq_1, fastq_2, out_file, outdir, ids, single_node=False
):
    if single_node:
        return align_sample_no_split(config, fastq_1, fastq_2, out_file, outdir, ids)
    else:
        return align_sample_split(config, fastq_1, fastq_2, out_file, outdir, ids)


def align_sample_no_split(config, fastq_1, fastq_2, out_file, outdir, ids):
    ref_genome = config['ref_genome']['file']

    markdups_metrics = os.path.join(outdir, 'markdups_metrics.pdf')
    samtools_flagstat = os.path.join(outdir, 'samtools_flagstat.txt')

    out_bai = out_file + '.bai'

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='align_bwa_mem',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00',
            ncpus=config['threads'],
            disk=300
        ),
        func=tasks.align_bwa_mem,
        args=(
            pypeliner.managed.InputFile(fastq_1),
            pypeliner.managed.InputFile(fastq_2),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam'),
            config['threads'],
            mgd.TempSpace("temp_align")
        ),
        kwargs={
            'sample_id': ids[0],
            'lane_id': ids[1],
            'read_group_info': config['read_group_info'],
            'docker_config': config['docker']
        }
    )

    workflow.transform(
        name='sort',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00',
            ncpus=config['threads'],
            disk=300
        ),
        func='wgs.workflows.alignment.tasks.bam_sort',
        args=(
            pypeliner.managed.TempInputFile('aligned.bam'),
            pypeliner.managed.TempOutputFile('sorted.bam'),
        ),
        kwargs={
            'docker_image': config['docker']['samtools'],
            'threads': config['threads'],
            # 'mem': '8G',
        }
    )

    workflow.transform(
        name='markdups',
        ctx=helpers.get_default_ctx(
            memory=12,
            walltime='24:00',
            ncpus=1,
            disk=300
        ),
        func=tasks.markdups,
        args=(
            pypeliner.managed.TempInputFile('sorted.bam'),
            pypeliner.managed.OutputFile(out_file),
            pypeliner.managed.OutputFile(markdups_metrics),
            pypeliner.managed.TempSpace("temp_markdups"),
        ),
        kwargs={'docker_image': config['docker']['picard']}
    )

    workflow.transform(
        name='index_and_flagstat',
        func=tasks.index_and_flagstat,
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='24:00',
            disk=200
        ),
        args=(
            pypeliner.managed.InputFile(out_file),
            pypeliner.managed.OutputFile(out_bai),
            pypeliner.managed.OutputFile(samtools_flagstat)
        ),
        kwargs={'docker_image': config['docker']['samtools']}
    )

    return workflow


def align_sample_split(config, fastq_1, fastq_2, out_file, outdir, ids):
    ref_genome = config['ref_genome']['file']

    read_group_config = config.get('read_group', {})

    markdups_metrics = os.path.join(outdir, 'markdups_metrics.pdf')
    samtools_flagstat = os.path.join(outdir, 'samtools_flagstat.txt')

    out_bai = out_file + '.bai'

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=pypeliner.managed.TempOutputObj('read_group_config'),
        value=read_group_config,
    )

    workflow.transform(
        name='split_fastq_1',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='12:00',
        ),
        func=biowrappers.components.io.fastq.tasks.split_fastq,
        args=(
            pypeliner.managed.InputFile(fastq_1),
            pypeliner.managed.TempOutputFile('read_1', 'split'),
            config['split_size'],
        ),
    )

    workflow.transform(
        name='split_fastq_2',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='12:00',
        ),
        func=biowrappers.components.io.fastq.tasks.split_fastq,
        args=(
            pypeliner.managed.InputFile(fastq_2),
            pypeliner.managed.TempOutputFile('read_2', 'split', axes_origin=[]),
            config['split_size'],
        ),
    )

    workflow.transform(
        name='align_bwa_mem',
        axes=('split',),
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='8:00',
            ncpus=config['threads'],
        ),
        func=tasks.align_bwa_mem,
        args=(
            pypeliner.managed.TempInputFile('read_1', 'split'),
            pypeliner.managed.TempInputFile('read_2', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
            config['threads'],
            mgd.TempSpace("temp_align_split", 'split')
        ),
        kwargs={
            'sample_id': ids[0],
            'lane_id': ids[1],
            'read_group_info': config['read_group_info'],
            'docker_config': config['docker']
        }
    )

    workflow.transform(
        name='sort',
        axes=('split',),
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='8:00',
        ),
        func='wgs.workflows.alignment.tasks.bam_sort',
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
        ),
        kwargs={
            'docker_image': config['docker']['samtools'],
        }
    )

    workflow.transform(
        name='merge',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func="wgs.workflows.alignment.tasks.merge_bams",
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.TempOutputFile('merged.bam'),
        ),
        kwargs={
            'picard_docker_image': config['docker']['picard'],
            'samtools_docker_image': config['docker']['samtools']
        }
    )

    workflow.transform(
        name='markdups',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='24:00',
        ),
        func=tasks.markdups,
        args=(
            pypeliner.managed.TempInputFile('merged.bam'),
            pypeliner.managed.OutputFile(out_file),
            pypeliner.managed.OutputFile(markdups_metrics),
            pypeliner.managed.TempSpace("temp_markdups"),
        ),
        kwargs={
            'mem': '8G',
            'docker_image': config['docker']['picard']
        },
    )

    workflow.commandline(
        name='index',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='8:00',
            docker_image=config['docker']['samtools']
        ),
        args=(
            'samtools',
            'index',
            pypeliner.managed.InputFile(out_file),
            pypeliner.managed.OutputFile(out_bai)
        ),
    )

    workflow.commandline(
        name='flagstat',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='8:00',
            docker_image=config['docker']['samtools']
        ),
        args=(
            'samtools',
            'flagstat',
            pypeliner.managed.InputFile(out_file),
            '>',
            pypeliner.managed.OutputFile(samtools_flagstat)
        ),
    )

    return workflow
