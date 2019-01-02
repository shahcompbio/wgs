import pypeliner
import os

def markdups(input, output, metrics, tempdir):
    cmd = ['picard', '-Xmx4G', '-Xms4G',
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


def merge_bams(inputs, output, mem="2G", **kwargs):
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
    merge_bams(inputs, output, docker_image=containers['picard'])
    bam_index(output, output_index, docker_image=containers['samtools'])
