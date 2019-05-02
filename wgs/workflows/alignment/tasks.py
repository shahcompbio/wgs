import os

import pypeliner


def markdups(input, output, metrics, tempdir, mem="2G"):
    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MarkDuplicates',
           'INPUT=' + input,
           'OUTPUT=' + output,
           'METRICS_FILE=' + metrics,
           'REMOVE_DUPLICATES=False',
           'ASSUME_SORTED=True',
           'VALIDATION_STRINGENCY=LENIENT',
           'TMP_DIR=' + tempdir,
           'MAX_RECORDS_IN_RAM=150000'
           ]

    pypeliner.commandline.execute(*cmd)


def picard_merge_bams(inputs, output, mem="2G", **kwargs):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000'
           ]

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    pypeliner.commandline.execute(*cmd, **kwargs)


def bam_index(infile, outfile, **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        **kwargs)


def merge_bams(inputs, output, containers):
    output_index = output + '.bai'
    if not containers:
        containers = {}
    picard_merge_bams(inputs, output, docker_image=containers.get('picard'))
    bam_index(output, output_index, docker_image=containers.get('samtools'))


def align_bwa_mem(read_1, read_2, ref_genome, aligned_bam, threads, sample_id=None, lane_id=None, read_group_info=None):
    if read_group_info:
        rg_id = read_group_info['ID'].format(**{'sample_id': sample_id, 'lane_id': lane_id})
        read_group = ['@RG', 'ID:{0}'.format(rg_id)]
        for key, value in sorted(read_group_info.items()):
            if key == 'ID':
                continue
            read_group.append(':'.join((key, value)))
        read_group = '\\t'.join(read_group)
    elif sample_id or lane_id:
        sample_id = sample_id if sample_id else ''
        lane_id = lane_id if lane_id else ''
        ids = '-'.join([sample_id, lane_id])
        read_group = ['@RG', 'ID:{0}'.format(ids)]
        read_group = '\\t'.join(read_group)
    else:
        read_group = None

    align_cmd = ['bwa', 'mem']
    if read_group:
        align_cmd.extend(['-R', read_group])
    align_cmd.extend(['-M', ref_genome, read_1, read_2, '-t', threads])

    to_bam_cmd = ['samtools', 'view', '-bSh', '-', '>', aligned_bam]

    cmd = align_cmd + ['|'] + to_bam_cmd

    pypeliner.commandline.execute(*cmd)
