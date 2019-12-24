import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows.alignment.dtypes import dtypes


def collect_bam_metrics(
        config, bam, markdups_metrics, outdir, metrics, sample_id
):
    '''
    calculates bam metrics in bams
    1. picard insert metrics
    2. picard GC metrics
    3. picard wgs metrics
    4. fastqc metrics

    :param config: config containing docker 
    images for metrics
    :param bams: sample:bam dictionary
    :param metrics_csv: output csv containing
        metrics
    :param single_node:
    '''

    workflow = pypeliner.workflow.Workflow()

    picard_insert_metrics = os.path.join(outdir, 'picard_insert_metrics.txt')
    picard_insert_pdf = os.path.join(outdir, 'picard_insert.pdf')
    flagstat_metrics = os.path.join(outdir, 'flagstat_metrics.txt')
    picard_GC_metrics = os.path.join(outdir, 'picard_GC_metrics.txt')
    picard_GC_summary = os.path.join(outdir, 'picard_GC_summary.txt')
    picard_GC_pdf = os.path.join(outdir, 'picard_GC.pdf')
    picard_wgs_metrics = os.path.join(outdir, 'picard_wgs_metrics.txt')

    workflow.transform(
        name="calc_picard_insert_metrics",
        func='wgs.workflows.alignment.tasks.bam_collect_insert_metrics',
        args=(
            mgd.InputFile(bam),
            mgd.OutputFile(flagstat_metrics),
            mgd.OutputFile(picard_insert_metrics),
            mgd.OutputFile(picard_insert_pdf),
            mgd.TempSpace('picard_insert'),
        ),
        kwargs={
            'picard_docker': config["docker"]["picard"],
            'samtools_docker': config["docker"]["samtools"]
        }
    )

    workflow.transform(
        name="calc_picard_gc_metrics",
        func='wgs.workflows.alignment.tasks.bam_collect_gc_metrics',
        args=(
            mgd.InputFile(bam),
            config["ref_genome"]['file'],
            mgd.OutputFile(picard_GC_metrics),
            mgd.OutputFile(picard_GC_summary),
            mgd.OutputFile(picard_GC_pdf),
            mgd.TempSpace('picard_gc')
        ),
        kwargs={'docker_image': config["docker"]["picard"]}
    )

    workflow.transform(
        name="calc_picard_wgs_metrics",
        func='wgs.workflows.alignment.tasks.bam_collect_wgs_metrics',
        args=(
            mgd.InputFile(bam),
            config['ref_genome']['file'],
            mgd.OutputFile(picard_wgs_metrics),
            config['picard_wgs_params'],
            mgd.TempSpace('picard_wgs')
        ),
        kwargs={'docker_image': config["docker"]["picard"]}
    )

    workflow.transform(
        name='collect_metrics',
        func='wgs.workflows.alignment.tasks.bam_collect_all_metrics',
        args=(
            mgd.InputFile(flagstat_metrics),
            mgd.InputFile(picard_insert_metrics),
            mgd.InputFile(picard_wgs_metrics),
            mgd.InputFile(markdups_metrics),
            mgd.OutputFile(metrics),
            sample_id
        ),
        kwargs= {
            'main_dtypes':dtypes()['metrics'],
            'insert_dtypes': dtypes()['insert_metrics']
        }
    )

    return workflow


def fastqc_workflow(fastq_r1, fastq_r2, outdir, config):
    report_r1 = os.path.join(outdir, 'R1_fastqc_report')
    report_r2 = os.path.join(outdir, 'R2_fastqc_report')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="fastqc_r1",
        func='wgs.workflows.alignment.tasks.run_fastqc',
        args=(
            mgd.InputFile(fastq_r1),
            mgd.OutputFile(report_r1),
            mgd.TempSpace('fastqc_R1'),
        ),
        kwargs={
            'docker_image': config["docker"]["fastqc"],
        }
    )

    workflow.transform(
        name="fastqc_r2",
        func='wgs.workflows.alignment.tasks.run_fastqc',
        args=(
            mgd.InputFile(fastq_r2),
            mgd.OutputFile(report_r2),
            mgd.TempSpace('fastqc_R2'),
        ),
        kwargs={
            'docker_image': config["docker"]["fastqc"],
        }
    )

    return workflow


def align_samples(
        config,
        config_globals,
        fastqs_r1,
        fastqs_r2,
        bam_outputs,
        outdir,
        single_node=False
):
    lane_metrics_template = os.path.join(outdir, '{sample_id}', 'metrics', 'lane_metrics', '{lane_id}')

    if single_node:
        align_func = align_sample_no_split
    else:
        align_func = align_sample_split

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'lane_id'),
        value=fastqs_r1.keys(),
    )

    workflow.subworkflow(
        name='fastqc_workflow',
        func=fastqc_workflow,
        axes=('sample_id', 'lane_id'),
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r2),
            mgd.Template('fastqc', 'sample_id', 'lane_id', template=lane_metrics_template),
            config
        )
    )

    workflow.subworkflow(
        name='align_samples',
        func=align_func,
        axes=('sample_id', 'lane_id'),
        args=(
            config,
            mgd.InputFile('input.r1.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'sample_id', 'lane_id', fnames=fastqs_r2),
            mgd.TempOutputFile('aligned_lanes.bam', 'sample_id', 'lane_id'),
            mgd.Template(lane_metrics_template, 'sample_id', 'lane_id'),
            [mgd.InputInstance("sample_id"),
             mgd.InputInstance("lane_id")]
        )
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
            mgd.TempOutputFile('merged_lanes.bam', 'sample_id', extensions=['.bai']),
        ),
        kwargs={
            'picard_docker_image': config['docker']['picard'],
            'samtools_docker_image': config['docker']['samtools']
        }
    )

    metrics_outdir = os.path.join(outdir, '{sample_id}', 'metrics')
    markdups_outputs = os.path.join(metrics_outdir, 'markdups_metrics.txt')
    metrics_output = os.path.join(outdir, '{sample_id}', '{sample_id}_metrics.csv')

    workflow.transform(
        name='markdups',
        ctx=helpers.get_default_ctx(
            memory=12,
            walltime='24:00',
            ncpus=1,
            disk=300
        ),
        func='wgs.workflows.alignment.tasks.markdups',
        axes=('sample_id',),
        args=(
            mgd.TempInputFile('merged_lanes.bam', 'sample_id', extensions=['.bai']),
            mgd.OutputFile('markdups.bam', 'sample_id', fnames=bam_outputs, extensions=['.bai']),
            mgd.OutputFile('markdups_metrics', 'sample_id', template=markdups_outputs),
            pypeliner.managed.TempSpace("temp_markdups", "sample_id"),
        ),
        kwargs={
            'picard_docker': config['docker']['picard'],
            'samtools_docker': config['docker']['samtools'],
        }
    )

    workflow.subworkflow(
        name='metrics',
        func=collect_bam_metrics,
        axes=('sample_id',),
        args=(
            config,
            mgd.InputFile('markdups.bam', 'sample_id', fnames=bam_outputs, extensions=['.bai']),
            mgd.InputFile('markdups_metrics', 'sample_id', template=markdups_outputs),
            mgd.Template('metrics_outdir', 'sample_id', template=metrics_outdir),
            mgd.OutputFile('metrics_output', 'sample_id', template=metrics_output),
            mgd.InputInstance('sample_id')
        )
    )

    return workflow


def align_sample_no_split(config, fastq_1, fastq_2, out_file, outdir, ids):
    ref_genome = config['ref_genome']['file']

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
        func='wgs.workflows.alignment.tasks.align_bwa_mem',
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
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'docker_image': config['docker']['samtools'],
            'threads': config['threads'],
            # 'mem': '8G',
        }
    )

    workflow.transform(
        name='index_and_flagstat',
        func='wgs.workflows.alignment.tasks.index_and_flagstat',
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
        func='biowrappers.components.io.fastq.tasks.split_fastq',
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
        func='biowrappers.components.io.fastq.tasks.split_fastq',
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
        func='wgs.workflows.alignment.tasks.align_bwa_mem',
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
            pypeliner.managed.OutputFile(out_file),
        ),
        kwargs={
            'picard_docker_image': config['docker']['picard'],
            'samtools_docker_image': config['docker']['samtools']
        }
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
