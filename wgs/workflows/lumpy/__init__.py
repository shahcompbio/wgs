import pypeliner
import pypeliner.managed as mgd

import tasks

def lumpy_preprocess_workflow(
        global_config, bamfile, sv_config, discordants_sorted_bam,
        splitters_sorted_bam, single_node=False
):

    workflow = pypeliner.workflow.Workflow()

    if single_node:
        workflow.transform(
            name='run_lumpy_preprocess',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpy_preprocess,
            args=(
                mgd.InputFile(bamfile),
                mgd.OutputFile(discordants_sorted_bam),
                mgd.OutputFile(splitters_sorted_bam),
                mgd.TempSpace("lumpy_preprocess_temp"),
                sv_config
            ),
        )
    else:
        workflow.transform(
            name='run_samtools_view_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_view,
            args=(
                mgd.InputFile(bamfile),
                mgd.TempOutputFile('normal.discordants.unsorted.bam'),
            ),
        )

        workflow.transform(
            name='run_lumpy_extract_split_reads_bwamem_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpy_extract_split_reads_bwamem,
            args=(
                mgd.InputFile(bamfile),
                mgd.TempOutputFile('normal.splitters.unsorted.bam'),
                sv_config,
            ),
        )

        workflow.transform(
            name='run_samtools_sort_discordants_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_sort,
            args=(
                mgd.TempInputFile('normal.discordants.unsorted.bam'),
                mgd.OutputFile(discordants_sorted_bam),
            ),
        )

        workflow.transform(
            name='run_samtools_sort_splitters_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_sort,
            args=(
                mgd.TempInputFile('normal.splitters.unsorted.bam'),
                mgd.OutputFile(splitters_sorted_bam),
            ),
        )

    return workflow


def create_lumpy_workflow(lumpy_vcf, global_config, sv_config, tumour_bam=None, normal_bam=None, single_node=False):
    workflow = pypeliner.workflow.Workflow()

    lumpy_job_name = 'run_lumpy'
    if normal_bam:
        normal_bam = mgd.InputFile(normal_bam)
        normal_disc = mgd.TempInputFile('normal.discordants.sorted.bam')
        normal_split = mgd.TempInputFile('normal.splitters.sorted.bam')
        lumpy_job_name += '_normal'
    else:
        normal_disc=None
        normal_split = None

    if tumour_bam:
        tumour_bam = mgd.InputFile(tumour_bam)
        tumour_disc = mgd.TempInputFile('tumour.discordants.sorted.bam')
        tumour_split = mgd.TempInputFile('tumour.splitters.sorted.bam')
        lumpy_job_name += '_tumour'
    else:
        tumour_disc=None
        tumour_split = None

    if normal_bam:
        workflow.subworkflow(
            name='preprocess_lumpy_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=lumpy_preprocess_workflow,
            args=(
                global_config,
                normal_bam,
                sv_config,
                mgd.TempOutputFile('normal.discordants.sorted.bam'),
                mgd.TempOutputFile('normal.splitters.sorted.bam')
            ),
            kwargs={'single_node': single_node}
        )

    if tumour_bam:
        workflow.subworkflow(
            name='preprocess_lumpy_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=lumpy_preprocess_workflow,
            args=(
                global_config,
                tumour_bam,
                sv_config,
                mgd.TempOutputFile('tumour.discordants.sorted.bam'),
                mgd.TempOutputFile('tumour.splitters.sorted.bam')
            ),
            kwargs={'single_node': single_node}
        )

    workflow.transform(
        name=lumpy_job_name,
        ctx={'mem': global_config['memory']['med'],
             'ncpus': 1},
        func=tasks.run_lumpyexpress,
        args=(
            mgd.OutputFile(lumpy_vcf),
            sv_config,
        ),
        kwargs={
            'tumour_bam': tumour_bam,
            'tumour_discordants': tumour_disc,
            'tumour_splitters': tumour_split,
            'normal_bam': normal_bam,
            'normal_discordants': normal_disc,
            'normal_splitters': normal_split
        }
    )

    return workflow
