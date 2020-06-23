'''
Created on Nov 1, 2015

@author: Andrew Roth
'''
from __future__ import division

import csv
import gzip
import os
import shutil

import numpy as np
import pypeliner
import pysam
from wgs.utils import helpers


def get_sample_id(bamfile):
    bam = pysam.AlignmentFile(bamfile)
    readgroups = bam.header['RG']

    samples = set()

    for readgroup in readgroups:
        samples.add(readgroup['SM'])

    assert len(samples) == 1

    return list(samples)[0]


def update_header_sample_ids(infile, outfile, tumour_id, normal_id):
    in_opener = gzip.open if '.gz' in infile else open
    out_opener = gzip.open if '.gz' in outfile else open

    with in_opener(infile) as indata:
        with out_opener(outfile, 'wt') as outdata:
            for line in indata:
                if line.startswith('#CHROM'):
                    outdata.write('##tumor_sample={}\n'.format(tumour_id))
                    outdata.write('##normal_sample={}\n'.format(normal_id))
                    line = line.replace('TUMOR', tumour_id).replace('NORMAL', normal_id)
                    outdata.write(line)
                else:
                    outdata.write(line)


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size) + 1)
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def get_chromosome_depth(chrom, bam_file, ref_genome, out_file, docker_image=None):
    chrom = chrom.split('_')[0]

    cmd = [
        'GetChromDepth',
        '--align-file', bam_file,
        '--chrom', chrom,
        '--output-file', out_file,
        # '--ref', ref_genome,
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def merge_chromosome_depths_weighted(infiles, outfile):
    data = {}

    for region, infile in infiles.items():
        size = int(region.split('_')[2]) - int(region.split('_')[1])
        with open(infile) as indata:
            depthdata = indata.readline()
            chrom, depth = depthdata.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append((float(depth) * size, size))

    with open(outfile, 'w') as output:
        for chrom, depths in data.items():
            total = sum([v[0] for v in depths])
            size = sum([v[1] for v in depths])
            output.write('{}\t{}\n'.format(chrom, total / size))


def merge_chromosome_depths_plain(infiles, outfile):
    data = {}

    if isinstance(infiles, dict):
        infiles = infiles.values()

    for infile in infiles:
        with open(infile) as indata:
            depthdata = indata.readline()
            chrom, depth = depthdata.strip().split()
            if chrom not in data:
                data[chrom] = []
            data[chrom].append(float(depth))

    with open(outfile, 'w') as output:
        for chrom, depths in data.items():
            output.write('{}\t{}\n'.format(chrom, np.mean(depths)))


def merge_chromosome_depths(infiles, outfile):
    if isinstance(infiles, dict):
        merge_chromosome_depths_weighted(infiles, outfile)
    else:
        merge_chromosome_depths_plain(infiles, outfile)


def strelka_one_node(
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_file,
        snv_file,
        tmp_dir,
        regions,
        known_sizes,
        is_exome=False,
        strelka_docker_image=None,
        vcftools_docker_image=None
):
    commands = []

    chromosomes = [val.split('_')[0] for val in regions]

    for chrom in chromosomes:
        chrom_temp_dir = os.path.join(tmp_dir, 'chroms', str(chrom))

        helpers.makedirs(chrom_temp_dir)

        outfile = os.path.join(chrom_temp_dir, 'depth.txt')

        cmd = [
            'GetChromDepth',
            '--align-file', normal_bam_file,
            '--chrom', chrom,
            '--output-file', outfile,
            # '--ref', ref_genome,
        ]

        commands.append(cmd)

    parallel_temp_dir = os.path.join(tmp_dir, 'gnu_parallel_temp_depths')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, strelka_docker_image)

    depthfiles = [os.path.join(tmp_dir, 'chroms', str(chrom), 'depth.txt') for chrom in chromosomes]
    depth_file = os.path.join(tmp_dir, 'chrom_depths.txt')
    merge_chromosome_depths_plain(depthfiles, depth_file)

    commands = []
    for i, region in enumerate(regions):
        ival_temp_dir = os.path.join(tmp_dir, 'intervals', str(i))
        helpers.makedirs(ival_temp_dir)
        indel_out = os.path.join(ival_temp_dir, 'strelka_indel.vcf')
        snv_out = os.path.join(ival_temp_dir, 'strelka_snv.vcf')
        stats_out = os.path.join(ival_temp_dir, 'stats.txt')

        cmd = genome_segment_cmd(
            depth_file,
            normal_bam_file,
            tumour_bam_file,
            ref_genome_fasta_file,
            indel_out,
            snv_out,
            stats_out,
            region,
            known_sizes,
            is_exome=is_exome,
        )
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tmp_dir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir, strelka_docker_image)

    indel_files = [os.path.join(tmp_dir, 'intervals', str(i), 'strelka_indel.vcf')
                   for i, region in enumerate(regions)]

    merge_temp = os.path.join(tmp_dir, 'snv_merge')
    snv_files = [os.path.join(tmp_dir, 'intervals', str(i), 'strelka_snv.vcf')
                 for i, region in enumerate(regions)]

    temp_strelka_snv = os.path.join(tmp_dir, 'snv_merge', 'temp_strelka_merge_snv.vcf')
    concatenate_vcf(snv_files, temp_strelka_snv, merge_temp, docker_image=vcftools_docker_image)

    temp_strelka_indel = os.path.join(tmp_dir, 'indel_merge' 'temp_strelka_merge_indel.vcf')
    concatenate_vcf(indel_files, temp_strelka_indel, merge_temp, docker_image=vcftools_docker_image)

    tumour_id = get_sample_id(tumour_bam_file)
    normal_id = get_sample_id(normal_bam_file)
    update_header_sample_ids(temp_strelka_snv, snv_file, tumour_id, normal_id)
    update_header_sample_ids(temp_strelka_indel, indel_file, tumour_id, normal_id)


def call_genome_segment(
        chrom_depth_file,
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_file,
        snv_file,
        tmp_dir,
        region,
        known_sizes,
        is_exome=False,
        docker_image=None,
):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    tmp_indel_file = os.path.join(tmp_dir, 'indels.vcf')

    tmp_snv_file = os.path.join(tmp_dir, 'snvs.vcf')

    stats_file = os.path.join(tmp_dir, 'stats.txt')

    cmd = genome_segment_cmd(
        chrom_depth_file,
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        tmp_indel_file,
        tmp_snv_file,
        stats_file,
        region,
        known_sizes,
        is_exome=is_exome,
    )

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    tumour_id = get_sample_id(tumour_bam_file)
    normal_id = get_sample_id(normal_bam_file)

    update_header_sample_ids(tmp_snv_file, snv_file, tumour_id, normal_id)
    update_header_sample_ids(tmp_indel_file, indel_file, tumour_id, normal_id)


def genome_segment_cmd(
        chrom_depth_file,
        normal_bam_file,
        tumour_bam_file,
        ref_genome_fasta_file,
        indel_file,
        snv_file,
        stats_file,
        region,
        known_sizes,
        is_exome=False,
        depthFilterMultiple=3.0,
        snvMaxFilteredBasecallFrac=0.4,
        snvMaxSpanningDeletionFrac=0.75,
        indelMaxWindowFilteredBasecallFrac=0.3,
        ssnvPrior=0.0001,
        sindelPrior=0.000001,
        ssnvNoise=0.0000000005,
        sindelNoiseFactor=2.2,
        ssnvNoiseStrandBiasFrac=0.0,
        minTier1Mapq=20,
        minTier2Mapq=0,
        ssnvQuality_LowerBound=15,
        sindelQuality_LowerBound=40,
        ssnvContamTolerance=0.15,
        indelContamTolerance=0.15,
):
    region = '{}:{}-{}'.format(*region.split('_'))

    genome_size = sum(known_sizes.values())

    cmd = [
        'run_strelka',
        normal_bam_file,
        tumour_bam_file,
        indel_file,
        snv_file,
        stats_file,
        # strelkaSharedWorkflow.py
        region,
        ref_genome_fasta_file,
        genome_size,
        '-max-indel-size', 50,
        # strelkaSomaticWorkflow.py
        '-min-mapping-quality', minTier1Mapq,
        '-min-qscore', 0,
        '-max-window-mismatch', 3, 20,
        '-indel-nonsite-match-prob', 0.5,
        '--somatic-snv-rate', ssnvPrior,
        '--shared-site-error-rate', ssnvNoise,
        '--shared-site-error-strand-bias-fraction', ssnvNoiseStrandBiasFrac,
        '--somatic-indel-rate', sindelPrior,
        '--shared-indel-error-factor', sindelNoiseFactor,
        '--tier2-min-mapping-quality', minTier2Mapq,
        '--tier2-mismatch-density-filter-count', 10,
        '--tier2-indel-nonsite-match-prob', 0.25,
        '--tier2-include-singleton',
        '--tier2-include-anomalous',
        '--strelka-snv-max-filtered-basecall-frac', snvMaxFilteredBasecallFrac,
        '--strelka-snv-max-spanning-deletion-frac', snvMaxSpanningDeletionFrac,
        '--strelka-snv-min-qss-ref', ssnvQuality_LowerBound,
        '--strelka-indel-max-window-filtered-basecall-frac', indelMaxWindowFilteredBasecallFrac,
        '--strelka-indel-min-qsi-ref', sindelQuality_LowerBound,
        '--ssnv-contam-tolerance', ssnvContamTolerance,
        '--indel-contam-tolerance', indelContamTolerance,
    ]

    if not is_exome:
        cmd.extend([
            '--strelka-chrom-depth-file', chrom_depth_file,
            '--strelka-max-depth-factor', depthFilterMultiple,
        ])

    return cmd


def get_known_chromosome_sizes(size_file, chromosomes):
    sizes = {}

    with open(size_file, 'r') as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            if row['chrom'] not in chromosomes:
                continue

            sizes[row['chrom']] = int(row['known_size'])

    return sizes


def count_fasta_bases(ref_genome_fasta_file, out_file, docker_image=None):
    cmd = [
        'countFastaBases',
        ref_genome_fasta_file,
        '>',
        out_file
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def index_bcf(in_file, docker_image=None):
    """ Index a VCF or BCF file with bcftools.
    :param in_file: Path of file to index.
    :param index_file: Path of index file.
    """
    pypeliner.commandline.execute('bcftools', 'index', in_file, docker_image=docker_image)


def index_vcf(vcf_file, docker_image=None):
    """ Create a tabix index for a VCF file
    :param vcf_file: Path of VCF to create index for. Should compressed by bgzip.
    :param index_file: Path of index file.
    This is meant to be used from pypeliner so it does some name mangling to add .tmp to the index file.
    """

    pypeliner.commandline.execute('tabix', '-f', '-p', 'vcf', vcf_file, docker_image=docker_image)


def concatenate_vcf(
        in_files, out_file, tempdir, docker_image=None,
        allow_overlap=False):
    """ Fast concatenation of VCF file using `bcftools`.
    :param in_files: dict with values being files to be concatenated. Files will be concatenated based on sorted order of keys.
    :param out_file: path where output file will be written in VCF format.
    """
    if isinstance(in_files, dict):
        in_files = in_files.values()

    helpers.makedirs(tempdir)

    merged_file = os.path.join(tempdir, 'merged.vcf')
    if allow_overlap:
        cmd = ['bcftools', 'concat', '-a', '-O', 'z', '-o', merged_file]
    else:
        cmd = ['bcftools', 'concat', '-O', 'z', '-o', merged_file]

    cmd += in_files

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    # sort merged vcf file
    cmd = ['bcftools', 'sort', '-O', 'z', '-o', out_file, merged_file]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    index_vcf(out_file, docker_image=docker_image)
    index_bcf(out_file, docker_image=docker_image)


def filter_vcf(raw_vcf, filtered_vcf, docker_image=None):
    cmd = [
        'bcftools',
        'view',
        '-O', 'z',
        '-f', '.,PASS',
        '-o', filtered_vcf,
        raw_vcf,
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    index_vcf(filtered_vcf, docker_image=docker_image)
    index_bcf(filtered_vcf, docker_image=docker_image)
