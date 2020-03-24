import pypeliner
import pypeliner.managed as mgd

from wgs.utils import helpers

from wgs.config import config


def lumpy_preprocess_workflow(
        bamfile, discordants_sorted_bam,
        splitters_sorted_bam, single_node=False
):
    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='run_lumpy_preprocess',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='96:00',
                disk=300
            ),
            func='wgs.workflows.lumpy.tasks.run_lumpy_preprocess',
            args=(
                mgd.InputFile(bamfile),
                mgd.OutputFile(discordants_sorted_bam),
                mgd.OutputFile(splitters_sorted_bam),
                mgd.TempSpace("lumpy_preprocess_temp"),
            ),
            kwargs={
                'lumpy_docker_image': config.containers('lumpy'),
                'samtools_docker_image': config.containers('samtools')
            }
        )
    else:
        workflow.transform(
            name='run_samtools_view_normal',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='24:00',
            ),
            func='wgs.workflows.lumpy.tasks.run_samtools_view',
            args=(
                mgd.InputFile(bamfile),
                mgd.TempOutputFile('normal.discordants.unsorted.bam'),
            ),
            kwargs={'docker_image': config.containers('samtools')}
        )

        workflow.transform(
            name='run_lumpy_extract_split_reads_bwamem_normal',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='24:00',
            ),
            func='wgs.workflows.lumpy.tasks.run_lumpy_extract_split_reads_bwamem',
            args=(
                mgd.InputFile(bamfile),
                mgd.TempOutputFile('normal.splitters.unsorted.bam'),
                config.default_params('breakpoint_calling')['lumpy_paths']
            ),
            kwargs={'docker_image': config.containers('lumpy')}
        )

        workflow.transform(
            name='run_samtools_sort_discordants_normal',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='24:00',
            ),
            func='wgs.workflows.lumpy.tasks.run_samtools_sort',
            args=(
                mgd.TempInputFile('normal.discordants.unsorted.bam'),
                mgd.OutputFile(discordants_sorted_bam),
            ),
            kwargs={'docker_image': config.containers('samtools')}
        )

        workflow.transform(
            name='run_samtools_sort_splitters_normal',
            ctx=helpers.get_default_ctx(
                memory=10,
                walltime='24:00',
            ),
            func='wgs.workflows.lumpy.tasks.run_samtools_sort',
            args=(
                mgd.TempInputFile('normal.splitters.unsorted.bam'),
                mgd.OutputFile(splitters_sorted_bam),
            ),
            kwargs={'docker_image': config.containers('samtools')}
        )

    return workflow


def create_lumpy_workflow(lumpy_vcf, tumour_bam=None, normal_bam=None, single_node=False):
    workflow = pypeliner.workflow.Workflow()

    lumpy_job_name = 'run_lumpy'
    if normal_bam:
        normal_bam = mgd.InputFile(normal_bam)
        normal_disc = mgd.TempInputFile('normal.discordants.sorted.bam')
        normal_split = mgd.TempInputFile('normal.splitters.sorted.bam')
        lumpy_job_name += '_normal'
    else:
        normal_disc = None
        normal_split = None

    if tumour_bam:
        tumour_bam = mgd.InputFile(tumour_bam)
        tumour_disc = mgd.TempInputFile('tumour.discordants.sorted.bam')
        tumour_split = mgd.TempInputFile('tumour.splitters.sorted.bam')
        lumpy_job_name += '_tumour'
    else:
        tumour_disc = None
        tumour_split = None

    if normal_bam:
        workflow.subworkflow(
            name='preprocess_lumpy_normal',
            func=lumpy_preprocess_workflow,
            args=(
                normal_bam,
                mgd.TempOutputFile('normal.discordants.sorted.bam'),
                mgd.TempOutputFile('normal.splitters.sorted.bam')
            ),
            kwargs={'single_node': single_node}
        )

    if tumour_bam:
        workflow.subworkflow(
            name='preprocess_lumpy_tumour',
            func=lumpy_preprocess_workflow,
            args=(
                tumour_bam,
                mgd.TempOutputFile('tumour.discordants.sorted.bam'),
                mgd.TempOutputFile('tumour.splitters.sorted.bam')
            ),
            kwargs={'single_node': single_node}
        )

    workflow.transform(
        name=lumpy_job_name,
        ctx=helpers.get_default_ctx(
            memory=10,
            disk=500,
            walltime='72:00'
        ),
        func='wgs.workflows.lumpy.tasks.run_lumpyexpress',
        args=(
            mgd.OutputFile(lumpy_vcf),
            config.default_params('breakpoint_calling')['lumpy_paths']
        ),
        kwargs={
            'tumour_bam': tumour_bam,
            'tumour_discordants': tumour_disc,
            'tumour_splitters': tumour_split,
            'normal_bam': normal_bam,
            'normal_discordants': normal_disc,
            'normal_splitters': normal_split,
            'docker_image': config.containers('lumpy')
        }
    )

    return workflow
