import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers
from wgs.workflows.alignment import align_pipelines
from wgs.workflows.alignment.dtypes import dtypes


def collect_bam_metrics(
        bam, sample_id, refdir,
        metrics, metrics_tar, bam_tdf, picard_mem=8
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

    ref_genome = config.refdir_data(refdir)['paths']['reference']

    picard_wgs_params = config.default_params('alignment')['picard_wgs_params']

    reftype = config.refdir_data(refdir)['params']['reference_type']

    workflow = pypeliner.workflow.Workflow()

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
            mgd.InputFile(bam, extensions=['.bai']),
            mgd.TempOutputFile('markdups_temp_rm.bam', extensions=['.bai']),
            mgd.TempOutputFile('markdups_metrics.txt'),
            pypeliner.managed.TempSpace("temp_markdups"),
        ),
        kwargs={
            'picard_docker': config.containers('picard'),
            'samtools_docker': config.containers('samtools'),
            'mem': '{}G'.format(picard_mem),
        }
    )

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
            mgd.TempOutputFile('flagstat_metrics.txt'),
            mgd.TempOutputFile('picard_insert_metrics.txt'),
            mgd.TempOutputFile('picard_insert.pdf'),
            mgd.TempSpace('picard_insert'),
        ),
        kwargs={
            'picard_docker': config.containers('picard'),
            'samtools_docker': config.containers('samtools'),
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
            mgd.TempOutputFile('picard_gc_metrics.txt'),
            mgd.TempOutputFile('picard_gc_summary.txt'),
            mgd.TempOutputFile('picard_gc.pdf'),
            mgd.TempSpace('picard_gc')
        ),
        kwargs={
            'docker_image': config.containers('picard'),
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
            mgd.TempOutputFile('picard_wgs_metrics.txt'),
            picard_wgs_params,
            mgd.TempSpace('picard_wgs')
        ),
        kwargs={
            'docker_image': config.containers('picard'),
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
        kwargs={'docker_image': config.containers('igvtools')}

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
            mgd.TempInputFile('flagstat_metrics.txt'),
            mgd.TempInputFile('picard_insert_metrics.txt'),
            mgd.TempInputFile('picard_wgs_metrics.txt'),
            mgd.TempInputFile('markdups_metrics.txt'),
            mgd.OutputFile(metrics, extensions=['.yaml']),
            sample_id
        ),
        kwargs={
            'main_dtypes': dtypes()['metrics'],
            'insert_dtypes': dtypes()['insert_metrics']
        }
    )

    workflow.transform(
        name='tar',
        func='wgs.utils.helpers.make_tar_from_files',
        axes=('sample_id',),
        args=(
            mgd.OutputFile(metrics_tar),
            [
                mgd.TempInputFile('picard_insert_metrics.txt'),
                mgd.TempInputFile('picard_insert.pdf'),
                mgd.TempInputFile('flagstat_metrics.txt'),
                mgd.TempInputFile('picard_gc_metrics.txt'),
                mgd.TempInputFile('picard_gc_summary.txt'),
                mgd.TempInputFile('picard_gc.pdf'),
                mgd.TempInputFile('picard_wgs_metrics.txt'),
                mgd.TempInputFile('markdups_metrics.txt'),
            ],
            mgd.TempSpace('wgs_metrics')
        )
    )

    return workflow


def fastqc_workflow(
        fastq_r1,
        fastq_r2,
        fastq_tar
):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('lane_id'),
        value=list(fastq_r1.keys()),
    )

    workflow.transform(
        name="fastqc_r1",
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='48:00',
            disk=400
        ),
        axes=('lane_id',),
        func='wgs.workflows.alignment.tasks.run_fastqc',
        args=(
            mgd.InputFile('fastq_r1.fastq.gz', 'lane_id', fnames=fastq_r1),
            mgd.TempOutputFile('r1_fastqc.html', 'lane_id'),
            mgd.TempOutputFile('r1_fastqc.pdf', 'lane_id'),
            mgd.TempSpace('fastqc_R1', 'lane_id'),
        ),
        kwargs={
            'docker_image': config.containers("fastqc"),
        }
    )

    workflow.transform(
        name="fastqc_r2",
        func='wgs.workflows.alignment.tasks.run_fastqc',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='48:00',
            disk=400
        ),
        axes=('lane_id',),
        args=(
            mgd.InputFile('fastq_r2.fastq.gz', 'lane_id', fnames=fastq_r2),
            mgd.TempOutputFile('r2_fastqc.html', 'lane_id'),
            mgd.TempOutputFile('r2_fastqc.pdf', 'lane_id'),
            mgd.TempSpace('fastqc_R2', 'lane_id'),
        ),
        kwargs={
            'docker_image': config.containers('fastqc'),
        }
    )

    workflow.transform(
        name='tar_fastq',
        func='wgs.utils.helpers.make_tar_from_files',
        args=(
            mgd.OutputFile(fastq_tar),
            [
                mgd.TempInputFile('r1_fastqc.html', 'lane_id'),
                mgd.TempInputFile('r1_fastqc.pdf', 'lane_id'),
                mgd.TempInputFile('r1_fastqc.html', 'lane_id'),
                mgd.TempInputFile('r1_fastqc.pdf', 'lane_id'),
            ],
            mgd.TempSpace('fastqc_metrics')
        )
    )

    return workflow


def align_samples(
        fastqs_r1,
        fastqs_r2,
        bam_outputs,
        metrics_outputs,
        metrics_tar,
        fastqc_tar,
        bam_tdf,
        sample_info,
        refdir,
        sample_id,
        single_node=False,
        picard_mem=8,
):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('lane_id'),
        value=list(fastqs_r1.keys()),
    )

    workflow.subworkflow(
        name='fastqc_workflow',
        func=fastqc_workflow,
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'lane_id', fnames=fastqs_r2),
            mgd.OutputFile(fastqc_tar),
        )
    )

    workflow.subworkflow(
        name='align_samples',
        func=align_pipelines.align_sample,
        axes=('sample_id',),
        args=(
            mgd.InputFile('input.r1.fastq.gz', 'lane_id', fnames=fastqs_r1),
            mgd.InputFile('input.r2.fastq.gz', 'lane_id', fnames=fastqs_r2),
            mgd.OutputFile(bam_outputs, extensions=['.bai']),
            refdir,
            single_node,
            sample_id,
            sample_info
        ),
        kwargs={'picard_mem': picard_mem}
    )

    workflow.subworkflow(
        name='metrics',
        func=collect_bam_metrics,
        args=(
            mgd.InputFile(bam_outputs, extensions=['.bai']),
            sample_id,
            refdir,
            mgd.OutputFile(metrics_outputs, extensions=['.yaml']),
            mgd.OutputFile(metrics_tar),
            mgd.OutputFile(bam_tdf),
        ),
        kwargs={'picard_mem': picard_mem}
    )

    return workflow
