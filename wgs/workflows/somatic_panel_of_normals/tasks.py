import os

import pypeliner
import pysam
from wgs.utils import helpers
from wgs.utils import vcfutils


def get_sample_id(bamfile):
    bam = pysam.AlignmentFile(bamfile)
    readgroups = bam.header['RG']

    samples = set()

    for readgroup in readgroups:
        samples.add(readgroup['SM'])

    assert len(samples) == 1

    return list(samples)[0]


def mutect_tumour_only_cmd(reference, interval, normal_bam, vcf_out, germline_resource):
    interval = interval.split('_')

    interval = '{}:{}-{}'.format(interval[0], interval[1], interval[2])

    normal_sample_id = get_sample_id(normal_bam)

    cmd = [
        'gatk', 'Mutect2', '-R', reference, '-I', normal_bam,
        '-tumor', normal_sample_id,
        '--germline-resource', germline_resource,
        '--normal', normal_sample_id, '-O', vcf_out, '--intervals', interval
    ]

    return cmd


def run_mutect(vcf, reference, interval, normal_bam, germline_resource):
    cmd = mutect_tumour_only_cmd(
        reference, interval, normal_bam, vcf, germline_resource
    )
    pypeliner.commandline.execute(*cmd)


def run_mutect_one_job(
        tempdir, vcf, reference, intervals, normal_bam, germline_resource
):
    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output_vcf = os.path.join(ival_temp_dir, 'mutect.vcf.gz')
        cmd = mutect_tumour_only_cmd(reference, interval, normal_bam, output_vcf, germline_resource)
        commands.append(cmd)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir)

    vcf_files = [os.path.join(tempdir, str(i), 'mutect.vcf.gz') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'mutect_merge')
    helpers.makedirs(merge_tempdir)
    merge_vcfs(vcf_files, vcf, merge_tempdir)


def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)


def mutect_panel_of_normals(sample_vcfs, pon_vcf, tempdir):
    helpers.makedirs(tempdir)
    args_file = os.path.join(tempdir, 'file_list.args')

    with open(args_file, 'wt') as args:
        for samp_vcf in sample_vcfs.values():
            args.write('{}\n'.format(samp_vcf))

    cmd = [
        'gatk', 'CreateSomaticPanelOfNormals',
        '-vcfs', args_file,
        '-O', pon_vcf
    ]
    pypeliner.commandline.execute(*cmd)
