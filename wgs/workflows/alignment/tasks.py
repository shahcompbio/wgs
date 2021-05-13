import logging
import os
import shutil

import pypeliner
from wgs.utils import helpers
from wgs.workflows.alignment.collect_metrics import CollectMetrics


def produce_fastqc_report(fastq_filename, output_html, output_plots, temp_dir,
                          **kwargs):
    temp_out_dir = os.path.join(temp_dir, 'out')
    temp_tmp_dir = os.path.join(temp_dir, 'tmp')
    helpers.makedirs(temp_out_dir)
    helpers.makedirs(temp_tmp_dir)

    pypeliner.commandline.execute(
        'fastqc',
        '--outdir=' + temp_out_dir,
        '--dir=' + temp_tmp_dir,
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

    output_basename = os.path.join(temp_out_dir, fastq_basename)

    shutil.move(output_basename + '_fastqc.zip', output_plots)
    shutil.move(output_basename + '_fastqc.html', output_html)


def run_fastqc(fastq1, html_file, plot_file, tempdir):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    # empty fastq files
    if os.stat(fastq1).st_size < 100:
        return

    if not os.path.getsize(fastq1) == 0:
        produce_fastqc_report(fastq1, html_file, plot_file, tempdir)
    else:
        logging.getLogger("single_cell.align.tasks").warn(
            "fastq file %s is empty, skipping fastqc" % fastq1)


def index_and_flagstat(bamfile, indexfile, flagstatfile):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd)

    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd)


def flagstat(bamfile, flagstatfile):
    cmd = ['samtools', 'flagstat', bamfile, '>', flagstatfile]
    pypeliner.commandline.execute(*cmd)


def index(bamfile, indexfile):
    cmd = ['samtools', 'index', bamfile, indexfile]
    pypeliner.commandline.execute(*cmd)


def cleanup_header(current_header, fixed_header):
    current_header = open(current_header, 'rt').readlines()

    fixed_data = []
    written = set()

    for line in current_header:
        line = line.strip()
        if '@RG' in line or '@HD' in line or '@SQ' in line:
            fixed_data.append(line)
        elif line.startswith('@PG'):
            assert 'CL' in line
            caller = line[line.index('CL:'):]
            caller = tuple(caller.split()[:2])
            if caller not in written:
                fixed_data.append(line)
                written.add(caller)
        else:
            raise Exception(line)

    with open(fixed_header, 'wt') as outwriter:
        for line in fixed_data:
            outwriter.write(line)


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
           'MAX_RECORDS_IN_RAM=150000',
           'QUIET=true'
           ]

    pypeliner.commandline.execute(*cmd)

    index(output, output + '.bai')


def picard_merge_bams(inputs, output, tempdir, mem="2G"):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000',
           'TMP_DIR=' + tempdir,
           'QUIET=true'
           ]

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    pypeliner.commandline.execute(*cmd)


def bam_index(infile, outfile, **kwargs):
    pypeliner.commandline.execute(
        'samtools', 'index',
        infile,
        outfile,
        **kwargs)


def merge_bams(inputs, output, tempdir, mem=2):
    output_index = output + '.bai'
    picard_merge_bams(
        inputs, output, tempdir,
        mem='{}G'.format(mem)
    )
    bam_index(output, output_index)


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
    read_group = ['@RG', 'ID:' + id_str.format(sample_id=sample_id, lane_id=lane_id)]

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
        sample_id=None, lane_id=None
):
    if lane_id in sample_info:
        sample_info = sample_info[lane_id]

    readgroup = get_readgroup(sample_info, sample_id, lane_id)

    bwa_mem_paired_end(
        read_1, read_2, aligned_bam, ref_genome,
        readgroup, threads
    )


def bam_sort(bam_filename, sorted_bam_filename, tempdir, threads=1, mem="2G"):
    helpers.makedirs(tempdir)

    prefix = os.path.join(tempdir, 'samtools_sort')

    pypeliner.commandline.execute(
        'samtools', 'sort', '-@', threads, '-m', mem,
        bam_filename,
        '-o',
        sorted_bam_filename,
        '-T',
        prefix)


# taken from single cell
def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            config, tempdir, mem="2G"):
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
        'QUIET=true',
    )


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           mem="2G"):
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
        'QUIET=true',
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
                               mem="2G"):
    bam_flagstat(
        bam_filename,
        flagstat_metrics_filename,
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
        'QUIET=true',
    )


def bam_collect_all_metrics(
        flagstat, insert, wgs, markdups_metrics, output, sample_id, main_dtypes=None, insert_dtypes=None
):
    collmet = CollectMetrics(
        wgs, insert, flagstat, markdups_metrics, output, sample_id, main_dtypes, insert_dtypes
    )
    collmet.main()


def get_igvtools_count(input_bam, counts_file, reference):
    counts_file_no_tmp = counts_file[:-4]

    cmd = ['igvtools', 'count', input_bam, counts_file_no_tmp, reference]

    pypeliner.commandline.execute(*cmd)

    os.rename(counts_file_no_tmp, counts_file)
