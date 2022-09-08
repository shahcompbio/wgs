import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows.alignment.dtypes import dtypes


def collect_bam_metrics(
        bam, markdups_metrics, sample_id, refdir,
        metrics, picard_insert_metrics, picard_insert_pdf,
        flagstat_metrics, picard_gc_metrics, picard_gc_summary,
        picard_gc_pdf, picard_wgs_metrics, bam_tdf, picard_mem=8
):
    '''
    calculates bam metrics in bams
    1. picard insert metrics
    2. picard GC metrics
    3. picard wgs metrics
    4. fastqc metrics

    :param config: config
    images for metrics
    :param bams: sample:bam dictionary
    :param metrics_csv: output csv containing
        metrics
    :param single_node:
    '''

    ref_genome = config.refdir_data(refdir)['paths']['reference']

    picard_wgs_params = config.default_params('alignment')['picard_wgs_params']

    reftype = config.refdir_data(refdir)['params']['reference_type']

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="calc_picard_insert_metrics",
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='72:00',
            disk=400
        ),
        func='wgs.workflows.alignment.tasks.bam_collect_insert_metrics',
        args=(
            mgd.InputFile(bam),
            mgd.OutputFile(flagstat_metrics),
            mgd.OutputFile(picard_insert_metrics),
            mgd.OutputFile(picard_insert_pdf),
            mgd.TempSpace('picard_insert'),
        ),
        kwargs={
            'mem': '{}G'.format(picard_mem)
        }
    )

    workflow.transform(
        name="calc_picard_gc_metrics",
        func='wgs.workflows.alignment.tasks.bam_collect_gc_metrics',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='72:00',
            disk=400
        ),
        args=(
            mgd.InputFile(bam),
            ref_genome,
            mgd.OutputFile(picard_gc_metrics),
            mgd.OutputFile(picard_gc_summary),
            mgd.OutputFile(picard_gc_pdf),
            mgd.TempSpace('picard_gc')
        ),
        kwargs={
            'mem': '{}G'.format(picard_mem)
        }
    )

    workflow.transform(
        name="calc_picard_wgs_metrics",
        func='wgs.workflows.alignment.tasks.bam_collect_wgs_metrics',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='72:00',
            disk=400
        ),
        args=(
            mgd.InputFile(bam),
            ref_genome,
            mgd.OutputFile(picard_wgs_metrics),
            picard_wgs_params,
            mgd.TempSpace('picard_wgs')
        ),
        kwargs={
            'mem': '{}G'.format(picard_mem)
        }
    )

    workflow.transform(
        name='igvtools_tdf',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='16:00',
        ),
        func='wgs.workflows.alignment.tasks.get_igvtools_count',
        args=(
            pypeliner.managed.InputFile(bam),
            pypeliner.managed.OutputFile(bam_tdf),
            reftype
        ),
    )

    workflow.transform(
        name='collect_metrics',
        func='wgs.workflows.alignment.tasks.bam_collect_all_metrics',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='4:00',
            disk=400
        ),
        args=(
            mgd.InputFile(flagstat_metrics),
            mgd.InputFile(picard_insert_metrics),
            mgd.InputFile(picard_wgs_metrics),
            mgd.InputFile(markdups_metrics),
            mgd.OutputFile(metrics, extensions=['.yaml']),
            sample_id
        ),
        kwargs={
            'main_dtypes': dtypes()['metrics'],
            'insert_dtypes': dtypes()['insert_metrics']
        }
    )

    return workflow


def fastqc_workflow(fastq_r1, fastq_r2, r1_html, r1_plot, r2_html, r2_plot):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="fastqc_r1",
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='48:00',
            disk=400
        ),
        func='wgs.workflows.alignment.tasks.run_fastqc',
        args=(
            mgd.InputFile(fastq_r1),
            mgd.OutputFile(r1_html),
            mgd.OutputFile(r1_plot),
            mgd.TempSpace('fastqc_R1'),
        ),
    )

    workflow.transform(
        name="fastqc_r2",
        func='wgs.workflows.alignment.tasks.run_fastqc',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='48:00',
            disk=400
        ),
        args=(
            mgd.InputFile(fastq_r2),
            mgd.OutputFile(r2_html),
            mgd.OutputFile(r2_plot),
            mgd.TempSpace('fastqc_R2'),
        ),
    )

    return workflow


def align_samples(
        fastqs_r1,
        fastqs_r2,
        bam_outputs,
        metrics_outputs,
        metrics_tar,
        bam_tdf,
        sample_info,
        refdir,
        single_node=False,
        picard_mem=8,
):
    if single_node:
        align_func = align_sample_no_split
    else:
        align_func = align_sample_split

    raise Exception(sample_info)

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.TempOutputObj('sampleinfo'),
        value=sample_info
    )

    workflow.setobj(
        obj=mgd.OutputChunks('lane_id'),
        value=list(fastqs_r1.keys()),
    )

    workflow.subworkflow(
        name='fastqc_workflow',
        func=fastqc_workflow,
        axes=('lane_id',),
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'lane_id', fnames=fastqs_r2),
            mgd.TempOutputFile('fastqc_R1.html', 'lane_id'),
            mgd.TempOutputFile('fastqc_R1.pdf', 'lane_id'),
            mgd.TempOutputFile('fastqc_R2.html', 'lane_id'),
            mgd.TempOutputFile('fastqc_R2.pdf', 'lane_id'),
        )
    )

    workflow.subworkflow(
        name='align_samples',
        func=align_func,
        axes=('lane_id',),
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'lane_id', fnames=fastqs_r2),
            mgd.TempOutputFile('aligned_lanes.bam', 'lane_id'),
            mgd.TempOutputFile('samtools_flagstat.txt', 'lane_id'),
            mgd.InputInstance("lane_id"),
            mgd.TempInputObj('sampleinfo'),
            refdir
        ),
        kwargs={'picard_mem': picard_mem}
    )

    workflow.transform(
        name='merge_tumour_lanes',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='24:00',
            disk=400
        ),
        func="wgs.workflows.alignment.tasks.merge_bams",
        args=(
            mgd.TempInputFile('aligned_lanes.bam', 'lane_id'),
            mgd.TempOutputFile('merged_lanes.bam', extensions=['.bai']),
            mgd.TempSpace('merge_tumour_lanes_tempdir')
        ),
        kwargs={
            'mem': picard_mem
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
        func='wgs.workflows.alignment.tasks.markdups',
        args=(
            mgd.TempInputFile('merged_lanes.bam', extensions=['.bai']),
            mgd.OutputFile(bam_outputs, extensions=['.bai']),
            mgd.TempOutputFile('markdups_metrics'),
            pypeliner.managed.TempSpace("temp_markdups"),
        ),
        kwargs={
            'mem': '{}G'.format(picard_mem),
        }
    )

    workflow.subworkflow(
        name='metrics',
        func=collect_bam_metrics,
        args=(
            mgd.InputFile(bam_outputs, extensions=['.bai']),
            mgd.TempInputFile('markdups_metrics'),
            sample_info['SM'],
            refdir,
            mgd.OutputFile(metrics_outputs, extensions=['.yaml']),
            mgd.TempOutputFile('picard_insert_metrics.txt'),
            mgd.TempOutputFile('picard_insert_metrics.pdf'),
            mgd.TempOutputFile('flagstat_metrics.txt'),
            mgd.TempOutputFile('picard_gc_metrics.txt'),
            mgd.TempOutputFile('picard_gc_summary.txt'),
            mgd.TempOutputFile('picard_gc.pdf'),
            mgd.TempOutputFile('picard_wgs_metrics.txt'),
            mgd.OutputFile(bam_tdf),
        )
    )

    workflow.transform(
        name='tar',
        func='wgs.utils.helpers.make_tar_from_files',
        args=(
            mgd.OutputFile(metrics_tar),
            [
                mgd.TempInputFile('picard_insert_metrics.txt'),
                mgd.TempInputFile('picard_insert_metrics.pdf'),
                mgd.TempInputFile('flagstat_metrics.txt'),
                mgd.TempInputFile('picard_gc_metrics.txt'),
                mgd.TempInputFile('picard_gc_summary.txt'),
                mgd.TempInputFile('picard_gc.pdf'),
                mgd.TempInputFile('picard_wgs_metrics.txt'),
                mgd.TempInputFile('markdups_metrics'),
                mgd.TempInputFile('fastqc_R1.html', 'lane_id'),
                mgd.TempInputFile('fastqc_R1.pdf', 'lane_id'),
                mgd.TempInputFile('fastqc_R2.html', 'lane_id'),
                mgd.TempInputFile('fastqc_R2.pdf', 'lane_id'),
            ],
            mgd.TempSpace('wgs_metrics')
        )
    )
    return workflow


def align_sample_no_split(
        fastq_1, fastq_2, out_file,
        samtools_flagstat,
        lane_id, sample_info, refdir,
        picard_mem=None
):
    ref_genome = config.refdir_data(refdir)['paths']['reference']

    out_bai = out_file + '.bai'

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='align_bwa_mem',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='48:00',
            ncpus='8',
            disk=300
        ),
        func='wgs.workflows.alignment.tasks.align_bwa_mem',
        args=(
            pypeliner.managed.InputFile(fastq_1),
            pypeliner.managed.InputFile(fastq_2),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam'),
            '8',
            sample_info,
        ),
        kwargs={
            'sample_id': sample_info['SM'],
            'lane_id': lane_id,
        }
    )

    workflow.transform(
        name='sort',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='48:00',
            ncpus='8',
            disk=300
        ),
        func='wgs.workflows.alignment.tasks.bam_sort',
        args=(
            pypeliner.managed.TempInputFile('aligned.bam'),
            pypeliner.managed.OutputFile(out_file),
            pypeliner.managed.TempSpace('bam_sort_tempdir')
        ),
        kwargs={
            'threads': '8',
            'mem': '{}G'.format(picard_mem)
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
    )

    return workflow


def align_sample_split(
        fastq_1, fastq_2, out_file,
        samtools_flagstat,
        lane_id, sample_info, refdir,
        picard_mem=2
):
    ref_genome = config.refdir_data(refdir)['paths']['reference']

    split_size = config.default_params('alignment')['split_size']

    out_bai = out_file + '.bai'

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='split_fastq_1',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='24:00',
        ),
        func='biowrappers.components.io.fastq.tasks.split_fastq',
        args=(
            pypeliner.managed.InputFile(fastq_1),
            pypeliner.managed.TempOutputFile('read_1', 'split'),
            split_size,
        ),
    )

    workflow.transform(
        name='split_fastq_2',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='24:00',
        ),
        func='biowrappers.components.io.fastq.tasks.split_fastq',
        args=(
            pypeliner.managed.InputFile(fastq_2),
            pypeliner.managed.TempOutputFile('read_2', 'split', axes_origin=[]),
            split_size,
        ),
    )

    workflow.transform(
        name='align_bwa_mem',
        axes=('split',),
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='16:00',
            ncpus=8,
        ),
        func='wgs.workflows.alignment.tasks.align_bwa_mem',
        args=(
            pypeliner.managed.TempInputFile('read_1', 'split'),
            pypeliner.managed.TempInputFile('read_2', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
            '8',
            sample_info,
        ),
        kwargs={
            'sample_id': sample_info['SM'],
            'lane_id': lane_id,
        }
    )

    workflow.transform(
        name='sort',
        axes=('split',),
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='16:00',
        ),
        func='wgs.workflows.alignment.tasks.bam_sort',
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
            pypeliner.managed.TempSpace('bam_sort_by_split', 'split')
        ),
        kwargs={
            'mem': '{}G'.format(picard_mem)
        }
    )

    workflow.transform(
        name='merge',
        ctx=helpers.get_default_ctx(
            memory=8,
            walltime='72:00',
        ),
        func="wgs.workflows.alignment.tasks.merge_bams",
        args=(
            pypeliner.managed.TempInputFile('sorted.bam', 'split'),
            pypeliner.managed.OutputFile(out_file),
            pypeliner.managed.TempSpace('bam_merge_by_split')
        ),
        kwargs={
            'mem': picard_mem
        }
    )

    workflow.commandline(
        name='index',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='16:00',
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
            walltime='16:00',
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
