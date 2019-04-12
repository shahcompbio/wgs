import pypeliner
import pypeliner.managed as mgd

import tasks


def create_lumpy_workflow(lumpy_vcf, global_config, sv_config, tumour_bam=None, normal_bam=None):
    workflow = pypeliner.workflow.Workflow()

    if normal_bam:
        workflow.transform(
            name='run_samtools_view_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_view,
            args=(
                mgd.InputFile(normal_bam),
                mgd.TempOutputFile('normal.discordants.unsorted.bam'),
            ),
        )

        workflow.transform(
            name='run_lumpy_extract_split_reads_bwamem_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpy_extract_split_reads_bwamem,
            args=(
                mgd.InputFile(normal_bam),
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
                mgd.TempOutputFile('normal.discordants.sorted.bam'),
            ),
        )

        workflow.transform(
            name='run_samtools_sort_splitters_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_sort,
            args=(
                mgd.TempInputFile('normal.splitters.unsorted.bam'),
                mgd.TempOutputFile('normal.splitters.sorted.bam'),
            ),
        )

    if tumour_bam:
        workflow.transform(
            name='run_samtools_view_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_view,
            args=(
                mgd.InputFile(tumour_bam),
                mgd.TempOutputFile('tumour.discordants.unsorted.bam'),
            ),
        )

        workflow.transform(
            name='run_lumpy_extract_split_reads_bwamem_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpy_extract_split_reads_bwamem,
            args=(
                mgd.InputFile(tumour_bam),
                mgd.TempOutputFile('tumour.splitters.unsorted.bam'),
                sv_config,
            ),
        )

        workflow.transform(
            name='run_samtools_sort_discordants_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_sort,
            args=(
                mgd.TempInputFile('tumour.discordants.unsorted.bam'),
                mgd.TempOutputFile('tumour.discordants.sorted.bam'),
            ),
        )

        workflow.transform(
            name='run_samtools_sort_splitters_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_samtools_sort,
            args=(
                mgd.TempInputFile('tumour.splitters.unsorted.bam'),
                mgd.TempOutputFile('tumour.splitters.sorted.bam'),
            ),
        )

    if tumour_bam and not normal_bam:
        workflow.transform(
            name='run_lumpyexpress_unpaired_tumour',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpyexpress,
            args=(
                mgd.OutputFile(lumpy_vcf),
                sv_config,
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam),
                'tumour_discordants': mgd.TempInputFile('tumour.discordants.sorted.bam'),
                'tumour_splitters': mgd.TempInputFile('tumour.splitters.sorted.bam')
            }
        )
    elif not tumour_bam and normal_bam:
        workflow.transform(
            name='run_lumpyexpress_unpaired_normal',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpyexpress,
            args=(
                mgd.OutputFile(lumpy_vcf),
                sv_config,
            ),
            kwargs={
                'normal_bam': mgd.InputFile(normal_bam),
                'normal_discordants': mgd.TempInputFile('normal.discordants.sorted.bam'),
                'normal_splitters': mgd.TempInputFile('normal.splitters.sorted.bam')
            }
        )
    else:
        workflow.transform(
            name='run_lumpyexpress_paired',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1},
            func=tasks.run_lumpyexpress,
            args=(
                mgd.OutputFile(lumpy_vcf),
                sv_config,
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam),
                'tumour_discordants': mgd.TempInputFile('tumour.discordants.sorted.bam'),
                'tumour_splitters': mgd.TempInputFile('tumour.splitters.sorted.bam'),
                'normal_bam': mgd.InputFile(normal_bam),
                'normal_discordants': mgd.TempInputFile('normal.discordants.sorted.bam'),
                'normal_splitters': mgd.TempInputFile('normal.splitters.sorted.bam')
            }
        )

    return workflow
