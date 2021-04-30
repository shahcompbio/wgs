import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def align_sample(
        fastqs_r1, fastqs_r2, bam_outputs,
        refdir, single_node, sample_id, sample_info, picard_mem=8
):
    if single_node:
        align_func = align_sample_no_split
    else:
        align_func = align_sample_split

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('lane_id'),
        value=list(fastqs_r1.keys()),
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
            sample_id,
            mgd.InputInstance("lane_id"),
            sample_info,
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
            'picard_docker_image': config.containers('picard'),
            'samtools_docker_image': config.containers('samtools'),
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
            'picard_docker': config.containers('picard'),
            'samtools_docker': config.containers('samtools'),
            'mem': '{}G'.format(picard_mem),
        }
    )

    return workflow


def align_sample_no_split(
        fastq_1, fastq_2, out_file,
        samtools_flagstat, sample_id,
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
            'sample_id': sample_id,
            'lane_id': lane_id,
            'docker_image': config.containers('bwa')
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
            'docker_image': config.containers('picard'),
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
        kwargs={'docker_image': config.containers('samtools')}
    )

    return workflow


def align_sample_split(
        fastq_1, fastq_2, out_file,
        samtools_flagstat, sample_id,
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
            'sample_id': sample_id,
            'lane_id': lane_id,
            'docker_image': config.containers('bwa')
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
            'docker_image': config.containers('samtools'),
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
            'picard_docker_image': config.containers('picard'),
            'samtools_docker_image': config.containers('samtools'),
            'mem': picard_mem
        }
    )

    workflow.commandline(
        name='index',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='16:00',
            docker_image=config.containers('samtools')
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
            docker_image=config.containers('samtools')
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
