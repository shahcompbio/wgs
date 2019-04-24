import pypeliner.managed as mgd
import pysam
from pypeliner.workflow import Workflow

import strelkautils as utils
import tasks
import vcf_tasks
from strelkautils import default_chromosomes


def create_strelka_workflow(
        normal_bam,
        tumour_bam,
        ref_genome_fasta_file,
        indel_vcf_file,
        snv_vcf_file,
        global_config,
        varcall_config,
        chromosomes=default_chromosomes,
        use_depth_thresholds=True,
        single_node=False,
):
    workflow = Workflow()

    workflow.transform(
        name='generate_intervals',
        func=tasks.generate_intervals,
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '01:00'},
        ret=mgd.OutputChunks('interval'),
        args=(
            varcall_config['reference'],
            varcall_config['chromosomes']
        )
    )

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=varcall_config['chromosomes'],
    )

    workflow.transform(
        name='count_fasta_bases',
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '01:00'},
        func=tasks.count_fasta_bases,
        args=(
            ref_genome_fasta_file,
            mgd.TempOutputFile('ref_base_counts.tsv'),
        ),
    )

    workflow.transform(
        name='get_chrom_sizes',
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '01:00'},
        func=tasks.get_known_chromosome_sizes,
        ret=mgd.TempOutputObj('known_sizes'),
        args=(
            mgd.TempInputFile('ref_base_counts.tsv'),
        ),
    )

    if single_node:
        workflow.transform(
            name='call_somatic_variants',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': global_config['threads'], 'walltime': '08:00'},
            func=tasks.call_somatic_variants_one_job,
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
                mgd.TempInputObj('known_sizes'),
                ref_genome_fasta_file,
                mgd.InputChunks('interval'),
                chromosomes,
                mgd.TempSpace("strelka_single_node_run"),
            ),
        )
    else:
        workflow.transform(
            name='call_somatic_variants',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': global_config['threads'], 'walltime': '08:00'},
            axes=('interval',),
            func=tasks.call_somatic_variants,
            args=(
                mgd.InputFile(normal_bam, extensions=['.bai']),
                mgd.InputFile(tumour_bam, extensions=['.bai']),
                mgd.TempInputObj('known_sizes'),
                ref_genome_fasta_file,
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('somatic.indels.unfiltered.vcf.window', 'interval'),
                mgd.TempOutputFile('somatic.snvs.unfiltered.vcf', 'interval'),
                mgd.TempOutputFile('strelka.stats', 'interval'),
                mgd.InputInstance('interval'),
            ),
        )

        workflow.transform(
            name='add_indel_filters',
            axes=('chrom',),
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1, 'walltime': '01:00'},
            func=tasks.filter_indel_file_list,
            args=(
                mgd.TempInputFile('somatic.indels.unfiltered.vcf', 'interval', axes_origin=[]),
                mgd.TempInputFile('strelka.stats', 'interval', axes_origin=[]),
                mgd.TempInputFile('somatic.indels.unfiltered.vcf.window', 'interval', axes_origin=[]),
                mgd.TempOutputFile('somatic.indels.filtered.vcf', 'chrom'),
                mgd.InputInstance('chrom'),
                mgd.TempInputObj('known_sizes'),
                mgd.InputChunks('interval'),
            ),
            kwargs={'use_depth_filter': use_depth_thresholds}
        )

        workflow.transform(
            name='add_snv_filters',
            axes=('chrom',),
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1, 'walltime': '01:00'},
            func=tasks.filter_snv_file_list,
            args=(
                mgd.TempInputFile('somatic.snvs.unfiltered.vcf', 'interval', axes_origin=[]),
                mgd.TempInputFile('strelka.stats', 'interval', axes_origin=[]),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf', 'chrom'),
                mgd.InputInstance('chrom'),
                mgd.TempInputObj('known_sizes'),
                mgd.InputChunks('interval'),
            ),
            kwargs={'use_depth_filter': use_depth_thresholds}
        )

        workflow.transform(
            name='merge_indels',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1, 'walltime': '01:00'},
            func=vcf_tasks.concatenate_vcf,
            args=(
                mgd.TempInputFile('somatic.indels.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.indels.filtered.vcf.gz'),
            ),
        )

        workflow.transform(
            name='merge_snvs',
            ctx={'mem': global_config['memory']['med'],
                 'ncpus': 1, 'walltime': '01:00'},
            func=vcf_tasks.concatenate_vcf,
            args=(
                mgd.TempInputFile('somatic.snvs.filtered.vcf', 'chrom'),
                mgd.TempOutputFile('somatic.snvs.filtered.vcf.gz'),
            ),
        )

    workflow.transform(
        name='filter_indels',
        ctx={'mem': global_config['memory']['med'],
             'ncpus': 1, 'walltime': '01:00'},
        func=vcf_tasks.filter_vcf,
        args=(
            mgd.TempInputFile('somatic.indels.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.indels.passed.vcf'),
        ),
    )

    workflow.transform(
        name='filter_snvs',
        ctx={'mem': global_config['memory']['med'],
             'ncpus': 1, 'walltime': '01:00'},
        func=vcf_tasks.filter_vcf,
        args=(
            mgd.TempInputFile('somatic.snvs.filtered.vcf.gz'),
            mgd.TempOutputFile('somatic.snvs.passed.vcf'),
        ),
    )

    workflow.transform(
        name='finalise_indels',
        ctx={
            'ncpus': 1, 'walltime': '01:00'},
        func=vcf_tasks.finalise_vcf,
        args=(
            mgd.TempInputFile('somatic.indels.passed.vcf'),
            mgd.OutputFile(indel_vcf_file),
        ),
    )

    workflow.transform(
        name='finalise_snvs',
        ctx={
            'ncpus': 1, 'walltime': '01:00'},
        func=vcf_tasks.finalise_vcf,
        args=(
            mgd.TempInputFile('somatic.snvs.passed.vcf'),
            mgd.OutputFile(snv_vcf_file),
        ),
    )

    return workflow


def get_chromosomes(bam_file, chromosomes=None):
    chromosomes = _get_chromosomes(bam_file, chromosomes)

    return dict(zip(chromosomes, chromosomes))


def _get_chromosomes(bam_file, chromosomes=None):
    bam = pysam.Samfile(bam_file, 'rb')

    if chromosomes is None:
        chromosomes = bam.references

    else:
        chromosomes = chromosomes

    return [str(x) for x in chromosomes]


def get_coords(bam_file, chrom, split_size):
    coords = {}

    bam = pysam.Samfile(bam_file, 'rb')

    chrom_lengths = dict(zip(bam.references, bam.lengths))

    length = chrom_lengths[chrom]

    lside_interval = range(1, length + 1, split_size)

    rside_interval = range(split_size, length + split_size, split_size)

    for coord_index, (beg, end) in enumerate(zip(lside_interval, rside_interval)):
        coords[coord_index] = (beg, end)

    return coords
