import gzip
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


def update_header_sample_ids(infile, outfile, tumour_id, normal_id):
    in_opener = gzip.open if '.gz' in infile else open
    out_opener = gzip.open if '.gz' in outfile else open

    with in_opener(infile) as indata:
        with out_opener(outfile, 'wt') as outdata:
            for line in indata:
                if line.startswith('#CHROM'):
                    if tumour_id:
                        outdata.write('##tumor_sample={}\n'.format(tumour_id))
                        line = line.replace('TUMOUR', tumour_id)
                    if normal_id:
                        outdata.write('##normal_sample={}\n'.format(normal_id))
                        line = line.replace('NORMAL', normal_id)
                outdata.write(line)


def run_museq(
        out, log, reference, interval, museq_params, tempdir,
        tumour_bam=None, normal_bam=None, return_cmd=False,
        titan_mode=False
):
    '''
    Run museq script for all chromosomes and merge VCF files

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to the temporary output VCF file for the merged VCF files
    :param log: path to the log file
    :param config: path to the config YAML file
    '''

    helpers.makedirs(tempdir)
    tempout = os.path.join(tempdir, 'museq_output.vcf')

    if titan_mode:
        cmd = ['museq_het']
    else:
        cmd = ['museq']

    if tumour_bam:
        cmd.append('tumour:' + tumour_bam)
    if normal_bam:
        cmd.append('normal:' + normal_bam)

    interval = interval.split('_')
    if len(interval) == 1:
        interval = interval[0]
    elif len(interval) == 2:
        interval = interval[0] + ':' + interval[1]
    elif len(interval) == 3:
        interval = interval[0] + ':' + interval[1] + '-' + interval[2]
    else:
        raise Exception()

    cmd.extend(['reference:' + reference, '--out', tempout,
                '--log', log, '--interval', interval, '-v'])

    if not tumour_bam or not normal_bam:
        cmd.extend(['-s'])

    for key, val in museq_params.items():
        if isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)

    if return_cmd:
        return cmd
    else:
        pypeliner.commandline.execute(*cmd)

    tumour_id = get_sample_id(tumour_bam) if tumour_bam else None
    normal_id = get_sample_id(normal_bam) if normal_bam else None

    update_header_sample_ids(tempout, out, tumour_id, normal_id)


def run_museq_one_job(
        tempdir, museq_vcf, reference, intervals, museq_params,
        tumour_bam=None, normal_bam=None, titan_mode=False
):
    '''
    Run museq script for all chromosomes and merge VCF files

    :param tumour: path to tumour bam
    :param normal: path to normal bam
    :param out: path to the temporary output VCF file for the merged VCF files
    :param log: path to the log file
    :param config: path to the config YAML file
    '''

    commands = []
    for i, interval in enumerate(intervals):
        ival_temp_dir = os.path.join(tempdir, str(i))
        helpers.makedirs(ival_temp_dir)
        output = os.path.join(ival_temp_dir, 'museq.vcf')
        log = os.path.join(ival_temp_dir, 'museq.log')

        command = run_museq(
            output, log, reference, interval, museq_params, ival_temp_dir,
            tumour_bam=tumour_bam, normal_bam=normal_bam,
            return_cmd=True, titan_mode=titan_mode
        )

        commands.append(command)

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    helpers.run_in_gnu_parallel(commands, parallel_temp_dir)

    vcf_files = [os.path.join(tempdir, str(i), 'museq.vcf') for i in range(len(intervals))]
    merge_tempdir = os.path.join(tempdir, 'museq_merge')
    helpers.makedirs(merge_tempdir)
    temp_museq_vcf = os.path.join(merge_tempdir, 'temp_museq_merge.vcf')
    merge_vcfs(vcf_files, temp_museq_vcf, merge_tempdir)

    tumour_id = get_sample_id(tumour_bam)
    normal_id = get_sample_id(normal_bam)
    update_header_sample_ids(temp_museq_vcf, museq_vcf, tumour_id, normal_id)

def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)
