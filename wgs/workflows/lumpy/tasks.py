import os

import pypeliner
from wgs.utils import helpers

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')


def run_samtools_view(infile, outfile, docker_image=None):
    cmd = ['samtools', 'view', '-b', '-F', '1294', infile, '>', outfile]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_lumpy_extract_split_reads_bwamem(infile, outfile, docker_image=None):
    cmd = [
        'samtools', 'view', '-h', infile, '|',
        'extractSplitReads_BwaMem', '-i', 'stdin', '|',
        'samtools', 'view', '-Sb', '-',
        '>', outfile
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_samtools_sort(infile, outfile, docker_image=None):
    cmd = ['samtools', 'sort', infile, '-o', outfile]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_lumpyexpress(
        lumpy_vcf, config,
        normal_bam=None,
        tumour_bam=None,
        normal_discordants=None,
        tumour_discordants=None,
        normal_splitters=None,
        tumour_splitters=None,
        docker_image=None
):
    lumpyexpress = config['lumpyexpress']

    cmd = [lumpyexpress,
           '-B', ','.join([e for e in [normal_bam, tumour_bam] if e]),
           '-S', ','.join([e for e in [normal_splitters, tumour_splitters] if e]),
           '-D', ','.join([e for e in [normal_discordants, tumour_discordants] if e]),
           '-o', lumpy_vcf]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def run_lumpy_preprocess(bamfile, disc_reads, split_reads, tempdir, samtools_docker_image=None, lumpy_docker_image=None):
    helpers.makedirs(tempdir)

    # disc
    unsorted_disc = os.path.join(tempdir, 'discordants.unsorted.bam')
    run_samtools_view(bamfile, unsorted_disc, docker_image=samtools_docker_image)
    run_samtools_sort(unsorted_disc, disc_reads, docker_image=samtools_docker_image)

    unsorted_split = os.path.join(tempdir, 'splitters.unsorted.bam')
    run_lumpy_extract_split_reads_bwamem(bamfile, unsorted_split, docker_image=lumpy_docker_image)
    run_samtools_sort(unsorted_split, split_reads, docker_image=samtools_docker_image)
