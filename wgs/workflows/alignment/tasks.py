import logging
import os
import shutil

import collect_metrics
import pypeliner
from wgs.utils import helpers


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir,
                          **kwargs):
    helpers.makedirs(temp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_dir,
        fastq_filename,
        **kwargs)

    fastq_basename = os.path.basename(fastq_filename)
    if fastq_basename.endswith(".fastq.gz"):
        fastq_basename = fastq_basename[:-len(".fastq.gz")]
    elif fastq_basename.endswith(".fq.gz"):
        fastq_basename = fastq_basename[:-len(".fq.gz")]
    elif fastq_basename.endswith(".fq"):
        fastq_basename = fastq_basename[:-len(".fq")]
    elif fastq_basename.endswith(".fastq"):
        fastq_basename = fastq_basename[:-len(".fastq")]
    else:
        raise Exception("Unknown file type")

    output_basename = os.path.join(temp_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def run_fastqc(fastq1, html_file, plot_file, tempdir, docker_image=None):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    # empty fastq files
    if os.stat(fastq1).st_size < 100:
        return

    if not os.path.getsize(fastq1) == 0:
        produce_fastqc_report(fastq1, html_file, plot_file, tempdir,
                              docker_image=docker_image)
    else:
        logging.getLogger("single_cell.align.tasks").warn(
            "fastq file %s is empty, skipping fastqc" % fastq1)


def index_and_flagstat(bamfile, indexfile, flagstatfile, docker_image=None):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def flagstat(bamfile, flagstatfile, docker_image=None):
    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def index(bamfile, indexfile, docker_image=None):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def markdups(input, output, metrics, tempdir, mem="2G", picard_docker=None, samtools_docker=None):
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

    pypeliner.commandline.execute(*cmd, docker_image=picard_docker)

    index(output, output + '.bai', docker_image=samtools_docker)


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


def merge_bams(inputs, output, picard_docker_image=None, samtools_docker_image=None):
    output_index = output + '.bai'
    picard_merge_bams(inputs, output, docker_image=picard_docker_image)
    bam_index(output, output_index, docker_image=samtools_docker_image)


def bwa_mem_paired_end(fastq1, fastq2, output,
                       reference, readgroup,
                       numthreads,
                       **kwargs):
    """
    run bwa aln on both fastq files,
    bwa sampe to align, and convert to bam with samtools view
    """

    if not numthreads:
        numthreads = 1

    if not readgroup:
        pypeliner.commandline.execute(
            'bwa', 'mem', '-M',
            '-t', numthreads,
            reference, fastq1, fastq2,
            '|', 'samtools', 'view', '-bSh', '-',
            '>', output,
            **kwargs)
    else:
        try:
            readgroup_literal = '"' + readgroup + '"'
            pypeliner.commandline.execute(
                'bwa', 'mem', '-M', '-R', readgroup_literal,
                '-t', numthreads,
                reference, fastq1, fastq2,
                '|', 'samtools', 'view', '-bSh', '-',
                '>', output,
                **kwargs)
        except pypeliner.commandline.CommandLineException:
            pypeliner.commandline.execute(
                'bwa', 'mem', '-M', '-R', readgroup,
                '-t', numthreads,
                reference, fastq1, fastq2,
                '|', 'samtools', 'view', '-bSh', '-',
                '>', output,
                **kwargs)


def get_readgroup(sample_info, sample_id, lane_id):
    if not sample_id or not lane_id:
        raise Exception('sample and lane ids are required')

    id_str = sample_info.pop('ID') if sample_info else 'ID:{0}_{1}'
    read_group = ['@RG', 'ID:'+id_str.format(sample_id = sample_id, lane_id=lane_id)]

    if sample_info:
        for key, value in sorted(sample_info.items()):
            value = value.format(sample_id=sample_id, lane_id=lane_id)
            read_group.append(':'.join((key, value)))

    read_group = '\\t'.join(read_group)

    return read_group


def samtools_sam_to_bam(samfile, bamfile,
                        **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'view', '-bSh', samfile,
        '>', bamfile,
        **kwargs)


def align_bwa_mem(
        read_1, read_2, ref_genome, aligned_bam, threads, sample_info,
        sample_id=None, lane_id=None, docker_config=None
):

    if lane_id in sample_info:
        sample_info = sample_info[lane_id]

    readgroup = get_readgroup(sample_info, sample_id, lane_id)

    bwa_mem_paired_end(
        read_1, read_2, aligned_bam, ref_genome,
        readgroup, threads, docker_image=docker_config['bwa']
    )


def bam_sort(bam_filename, sorted_bam_filename, threads=1, mem="2G", docker_image=None):
    pypeliner.commandline.execute(
        'samtools', 'sort', '-@', threads, '-m', mem,
        bam_filename,
        '-o',
        sorted_bam_filename,
        docker_image=docker_image)


# taken from single cell
def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            config, tempdir, mem="2G", docker_image=None):
    helpers.makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectWgsMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'MINIMUM_BASE_QUALITY=' +
                  str(config['min_bqual']),
                  'MINIMUM_MAPPING_QUALITY=' +
                  str(config['min_mqual']),
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
                  'COUNT_UNPAIRED=' +
                  ('True' if config['count_unpaired'] else 'False'),
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           mem="2G", docker_image=None):
    helpers.makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectGcBiasMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'S=' + summary_filename,
                  'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=docker_image
    )


def bam_flagstat(bam, metrics, **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'flagstat',
        bam,
        '>',
        metrics,
        **kwargs)


def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename,
                               metrics_filename, histogram_filename, tempdir,
                               mem="2G", picard_docker=None, samtools_docker=None):
    bam_flagstat(
        bam_filename,
        flagstat_metrics_filename,
        docker_image=samtools_docker
    )

    # Check if any paired reads exist
    has_paired = None
    with open(flagstat_metrics_filename) as f:
        for line in f:
            if 'properly paired' in line:
                if line.startswith('0 '):
                    has_paired = False
                else:
                    has_paired = True

    if has_paired is None:
        raise Exception('Unable to determine number of properly paired reads from {}'.format(
            flagstat_metrics_filename))

    if not has_paired:
        with open(metrics_filename, 'w') as f:
            f.write('## FAILED: No properly paired reads\n')
        with open(histogram_filename, 'w'):
            pass
        return

    helpers.makedirs(tempdir)

    pypeliner.commandline.execute(
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectInsertSizeMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        docker_image=picard_docker
    )


def bam_collect_all_metrics(
        flagstat, insert, wgs, markdups_metrics, output, sample_id, main_dtypes=None, insert_dtypes=None
):
    collmet = collect_metrics.CollectMetrics(
        wgs, insert, flagstat, markdups_metrics, output, sample_id, main_dtypes, insert_dtypes
    )
    collmet.main()
