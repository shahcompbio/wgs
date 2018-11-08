'''
Created on Nov 1, 2015

@author: Andrew Roth
'''
from __future__ import division

import os
import csv
import math
import pandas as pd
import pypeliner
import re
import vcf
import pysam

FILTER_ID_BASE = 'BCNoise'
FILTER_ID_DEPTH = 'DP'
FILTER_ID_INDEL_HPOL = 'iHpol'
FILTER_ID_QSI = 'QSI_ref'
FILTER_ID_QSS = 'QSS_ref'
FILTER_ID_REPEAT = 'Repeat'
FILTER_ID_SPANNING_DELETION = 'SpanDel'

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')


def generate_intervals(ref, chromosomes, size=100000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue

        for i in range(int((length/size)+1)):
            intervals.append( name+ "_" + str(i*size) +"_"+ str((i+1)*size))

    return intervals


def _get_files_for_chrom(infiles, intervals, chrom):

    if not isinstance(infiles, dict):
        infiles = {ival:infiles[ival] for ival in intervals}

    outfiles = {}

    for interval in intervals:
        ival_chrom = interval.split("_")[0]

        if ival_chrom==chrom:
            outfiles[interval] = infiles[interval]

    return outfiles


def get_known_chromosome_sizes(size_file, chromosomes):
    sizes = {}

    with open(size_file, 'r') as fh:
        reader = csv.DictReader(fh, ['path', 'chrom', 'known_size', 'size'], delimiter='\t')

        for row in reader:
            if row['chrom'] not in chromosomes:
                continue

            sizes[row['chrom']] = int(row['known_size'])

    return sizes

def count_fasta_bases(ref_genome_fasta_file, out_file):

    cmd = [
        'countFastaBases',
        ref_genome_fasta_file,
        '>',
        out_file
    ]

    pypeliner.commandline.execute(*cmd)


def call_somatic_variants(
        normal_bam_file,
        tumour_bam_file,
        known_sizes,
        ref_genome,
        indel_file,
        indel_window_file,
        snv_file,
        stats_file,
        interval,
        max_input_depth=10000,
        min_tier_one_mapq=20,
        min_tier_two_mapq=5,
        sindel_noise=0.000001,
        sindel_prior=0.000001,
        ssnv_noise=0.0000005,
        ssnv_noise_strand_bias_frac=0.5,
        ssnv_prior=0.000001):

    chrom, beg, end = interval.split('_')

    known_chrom_sizes = known_sizes[chrom]
    beg = int(beg)
    beg = beg+1 if beg == 0 else beg
    end = int(end)

    cmd = [
        'strelka2',
        '-bam-file', normal_bam_file,
        '--tumor-bam-file', tumour_bam_file,
        '-samtools-reference', ref_genome,
        '--somatic-indel-file', indel_file,
        '--somatic-snv-file', snv_file,
        '--variant-window-flank-file', 50, indel_window_file,
        '-bam-seq-name', chrom,
        '-report-range-begin', beg,
        '-report-range-end', end,
        '-clobber',
        '-filter-unanchored',
        '-genome-size', known_chrom_sizes,
        '-indel-nonsite-match-prob', 0.5,
        '-max-indel-size', 50,
        '-max-window-mismatch', 3, 20,
        '-min-paired-align-score', min_tier_one_mapq,
        '-min-single-align-score', 10,
        '-min-qscore', 0,
        '-print-used-allele-counts',
        '--max-input-depth', max_input_depth,
        '--min-contig-open-end-support', 35,
        '--report-file', stats_file,
        '--shared-site-error-rate', ssnv_noise,
        '--shared-site-error-strand-bias-fraction', ssnv_noise_strand_bias_frac,
        '--somatic-indel-rate', sindel_prior,
        '--shared-indel-error-rate', sindel_noise,
        '--somatic-snv-rate', ssnv_prior,
        '--tier2-include-anomalous',
        '--tier2-include-singleton',
        '--tier2-indel-nonsite-match-prob', 0.25,
        '--tier2-min-paired-align-score', min_tier_two_mapq,
        '--tier2-min-single-align-score', min_tier_two_mapq,
        '--tier2-mismatch-density-filter-count', 10,
        '--tier2-no-filter-unanchored',
        '--tier2-single-align-score-rescue-mode'
    ]

    pypeliner.commandline.execute(*cmd)


#=======================================================================================================================
# SNV filtering
#=======================================================================================================================


def filter_snv_file_list(
        in_files,
        stats_files,
        out_file,
        chrom,
        known_chrom_size,
        intervals,
        depth_filter_multiple=3.0,
        max_filtered_basecall_frac=0.4,
        max_spanning_deletion_frac=0.75,
        quality_lower_bound=15,
        use_depth_filter=True):


    known_chrom_size = known_chrom_size[chrom]

    in_files = _get_files_for_chrom(in_files, intervals, chrom)
    stats_files = _get_files_for_chrom(stats_files, intervals, chrom)

    max_normal_coverage = _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files)

    writer = None

    with open(out_file, 'wb') as out_fh:
        for key in sorted(in_files):
            reader = vcf.Reader(filename=in_files[key])

            if writer is None:
                # Add filters to header
                if use_depth_filter:
                    reader.filters[FILTER_ID_DEPTH] = vcf.parser._Filter(
                        id=FILTER_ID_DEPTH,
                        desc='Greater than {0}x chromosomal mean depth in Normal sample'.format(depth_filter_multiple)
                    )

                reader.filters[FILTER_ID_BASE] = vcf.parser._Filter(
                    id=FILTER_ID_BASE,
                    desc='Fraction of basecalls filtered at this site in either sample is at or above {0}'.format(
                        max_filtered_basecall_frac)
                )

                reader.filters[FILTER_ID_SPANNING_DELETION] = vcf.parser._Filter(
                    id=FILTER_ID_SPANNING_DELETION,
                    desc='Fraction of reads crossing site with spanning deletions in either sample exceeeds {0}'.format(
                        max_spanning_deletion_frac)
                )

                reader.filters[FILTER_ID_QSS] = vcf.parser._Filter(
                    id=FILTER_ID_QSS,
                    desc='Normal sample is not homozygous ref or ssnv Q-score < {0}, ie calls with NT!=ref or QSS_NT < {0}'.format(
                        quality_lower_bound)
                )

                writer = vcf.Writer(out_fh, reader)

            for record in reader:
                normal = record.genotype('NORMAL')

                tumour = record.genotype('TUMOR')

                # Normal depth filter
                if use_depth_filter and (normal.data.DP > max_normal_coverage):
                    record.add_filter(FILTER_ID_DEPTH)

                # Filtered basecall fraction
                normal_filtered_base_call_fraction = _get_filtered_base_call_fraction(normal.data)

                tumour_filtered_base_call_fraction = _get_filtered_base_call_fraction(tumour.data)

                if (normal_filtered_base_call_fraction >= max_filtered_basecall_frac) or \
                   (tumour_filtered_base_call_fraction >= max_filtered_basecall_frac):

                    record.add_filter(FILTER_ID_BASE)

                # Spanning deletion fraction
                normal_spanning_deletion_fraction = _get_spanning_deletion_fraction(normal.data)

                tumour_spanning_deletion_fraction = _get_spanning_deletion_fraction(tumour.data)

                if (normal_spanning_deletion_fraction > max_spanning_deletion_frac) or \
                   (tumour_spanning_deletion_fraction > max_spanning_deletion_frac):

                    record.add_filter(FILTER_ID_SPANNING_DELETION)

                # Q-val filter
                if (record.INFO['NT'] != 'ref') or (record.INFO['QSS_NT'] < quality_lower_bound):
                    record.add_filter(FILTER_ID_QSS)

                writer.write_record(record)

        if writer:
            writer.close()


def _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files):

    normal_coverage = _get_normal_coverage(stats_files.values())

    normal_mean_coverage = normal_coverage / known_chrom_size

    max_normal_coverage = normal_mean_coverage * depth_filter_multiple

    return max_normal_coverage


def _get_normal_coverage(stats_files):
    total_coverage = 0

    for file_name in stats_files:
        mean_matcher = re.compile('mean:\s(.*?)\s')

        sample_size_matcher = re.compile('sample_size:\s(.*?)\s')

        with open(file_name) as fh:
            for line in fh:
                if line.startswith('NORMAL_NO_REF_N_COVERAGE '):
                    mean = float(mean_matcher.search(line).group(1))

                    sample_size = float(sample_size_matcher.search(line).group(1))

                    if math.isnan(mean) or math.isnan(sample_size):
                        continue

                    total_coverage += mean * sample_size

    return total_coverage


def _get_filtered_base_call_fraction(data):
    frac = 0

    if data.DP > 0:
        frac = data.FDP / data.DP

    return frac


def _get_spanning_deletion_fraction(data):
    total = data.DP + data.SDP

    frac = 0

    if total > 0:
        frac = data.SDP / total

    return frac

#=======================================================================================================================
# Indel filtering
#=======================================================================================================================


def filter_indel_file_list(
        vcf_files,
        stats_files,
        window_files,
        out_file,
        chrom,
        known_chrom_size,
        intervals,
        depth_filter_multiple=3.0,
        max_int_hpol_length=14,
        max_ref_repeat=8,
        max_window_filtered_basecall_frac=0.3,
        quality_lower_bound=30,
        use_depth_filter=True):

    window_cols = (
        'chrom',
        'coord',
        'normal_window_used',
        'normal_window_filtered',
        'normal_window_submap',
        'tumour_window_used',
        'tumour_window_filtered',
        'tumour_window_submap'
    )


    known_chrom_size = known_chrom_size[chrom]

    vcf_files = _get_files_for_chrom(vcf_files, intervals, chrom)
    stats_files = _get_files_for_chrom(stats_files, intervals, chrom)
    window_files = _get_files_for_chrom(window_files, intervals, chrom)

    max_normal_coverage = _get_max_normal_coverage(chrom, depth_filter_multiple, known_chrom_size, stats_files)

    writer = None

    with open(out_file, 'wb') as out_fh:
        for key in sorted(vcf_files):
            window = pd.read_csv(
                window_files[key],
                comment='#',
                        converters={'chrom': str},
                header=None,
                names=window_cols,
                sep='\t')

            reader = vcf.Reader(filename=vcf_files[key])

            if writer is None:
                # Add format to header
                reader.formats['DP50'] = vcf.parser._Format(
                    id='DP50',
                    num=1,
                    type='Float',
                    desc='Average tier1 read depth within 50 bases'
                )

                reader.formats['FDP50'] = vcf.parser._Format(
                    id='FDP50',
                    num=1,
                    type='Float',
                    desc='Average tier1 number of basecalls filtered from original read depth within 50 bases'
                )

                reader.formats['SUBDP50'] = vcf.parser._Format(
                    id='SUBDP50',
                    num=1,
                    type='Float',
                    desc='Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases'
                )

                # Add filters to header
                if use_depth_filter:
                    reader.filters[FILTER_ID_DEPTH] = vcf.parser._Filter(
                        id=FILTER_ID_DEPTH,
                        desc='Greater than {0}x chromosomal mean depth in Normal sample'.format(depth_filter_multiple)
                    )

                reader.filters[FILTER_ID_REPEAT] = vcf.parser._Filter(
                    id=FILTER_ID_REPEAT,
                    desc='Sequence repeat of more than {0}x in the reference sequence'.format(max_ref_repeat)
                )

                reader.filters[FILTER_ID_INDEL_HPOL] = vcf.parser._Filter(
                    id=FILTER_ID_INDEL_HPOL,
                    desc='Indel overlaps an interrupted homopolymer longer than {0}x in the reference sequence'.format(
                        max_int_hpol_length)
                )

                reader.filters[FILTER_ID_BASE] = vcf.parser._Filter(
                    id=FILTER_ID_BASE,
                    desc='Average fraction of filtered basecalls within 50 bases of the indel exceeds {0}'.format(
                        max_window_filtered_basecall_frac)
                )

                reader.filters[FILTER_ID_QSI] = vcf.parser._Filter(
                    id=FILTER_ID_QSI,
                    desc='Normal sample is not homozygous ref or sindel Q-score < {0}, ie calls with NT!=ref or QSI_NT < {0}'.format(
                        quality_lower_bound)
                )

                writer = vcf.Writer(out_fh, reader)

            for record in reader:
                window_row = window.loc[
                    (window['chrom'] == str(record.CHROM)) & (window['coord'] == record.POS)].iloc[0]

                normal = record.genotype('NORMAL')

                tumour = record.genotype('TUMOR')

                normal_data = normal.data._asdict()

                tumour_data = tumour.data._asdict()

                # Add window data to vcf record
                record.add_format('DP50')

                record.add_format('FDP50')

                record.add_format('SUBDP50')

                normal_data['DP50'] = window_row['normal_window_used'] + window_row['normal_window_filtered']

                normal_data['FDP50'] = window_row['normal_window_filtered']

                normal_data['SUBDP50'] = window_row['normal_window_submap']

                tumour_data['DP50'] = window_row['tumour_window_used'] + window_row['tumour_window_filtered']

                tumour_data['FDP50'] = window_row['tumour_window_filtered']

                tumour_data['SUBDP50'] = window_row['tumour_window_submap']

                normal.data = _convert_dict_to_call(normal_data)

                tumour.data = _convert_dict_to_call(tumour_data)

                print normal.data

                # Add filters

                # Normal depth filter
                if use_depth_filter and (normal.data.DP > max_normal_coverage):
                    record.add_filter(FILTER_ID_DEPTH)

                # Ref repeat
                if 'RC' in record.INFO:
                    ref_repeat = record.INFO['RC']

                    if ref_repeat > max_ref_repeat:
                        record.add_filter(FILTER_ID_REPEAT)

                # Indel homopolymer
                if 'IHP' in record.INFO:
                    indel_homopolymer = record.INFO['IHP']

                    if indel_homopolymer > max_int_hpol_length:
                        record.add_filter(FILTER_ID_INDEL_HPOL)

                # Base filter
                normal_filtered_base_call_fraction = _get_filtered_base_call_fraction_indel(normal.data)

                tumour_filtered_base_call_fraction = _get_filtered_base_call_fraction_indel(tumour.data)

                if (normal_filtered_base_call_fraction >= max_window_filtered_basecall_frac) or \
                   (tumour_filtered_base_call_fraction >= max_window_filtered_basecall_frac):

                    record.add_filter(FILTER_ID_BASE)

                # Q-val filter
                if (record.INFO['NT'] != 'ref') or (record.INFO['QSI_NT'] < quality_lower_bound):
                    record.add_filter(FILTER_ID_QSI)

                writer.write_record(record)

        if writer:
            writer.close()


def _convert_dict_to_call(data_dict):
    call_data_class = vcf.model.make_calldata_tuple(data_dict.keys())

    return call_data_class(**data_dict)


def _get_filtered_base_call_fraction_indel(data):
    frac = 0

    if data.DP50 > 0:
        frac = data.FDP50 / data.DP50

    return frac


#=======================================================================================================================
# Indel and SNV Flagging
#=======================================================================================================================

def run_snpeff(infile, output, config):
    '''
    Run snpEff script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    :param config: path to the config YAML file
    '''

    snpeff_config = config['snpeff_params']['snpeff_config']
    genome_version = config['snpeff_params']['genome_version']

    cmd = ['snpEff', '-Xmx5G', '-Xms5G', '-XX:ParallelGCThreads=1',
           '-c', snpeff_config, genome_version, '-noStats', infile, '>',
           output]

    pypeliner.commandline.execute(*cmd)


def run_mutation_assessor(infile, output, config):
    '''
    Run Mutation Assessor script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    mutation_assessor_script = os.path.join(scripts_directory, 'annotate_mutation_assessor.py')
    db = config['mutation_assessor_params']['db']

    cmd = ['python', mutation_assessor_script, '--vcf', infile, '--output',
           output, '--db', db]

    pypeliner.commandline.execute(*cmd)


def run_DBSNP(infile, output, config):
    '''
    Run DBSNP script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'flagPos.py')
    db = config['dbsnp_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'DBSNP', '--input_type', 'snv', '--flag_with_id', 'True']

    pypeliner.commandline.execute(*cmd)


def run_1000gen(infile, output, config):
    '''
    Run 1000Gen script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'flagPos.py')
    db = config['1000gen_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', '1000Gen', '--input_type', 'snv']

    pypeliner.commandline.execute(*cmd)


def run_cosmic(infile, output, config):
    '''
    Run Cosmic script on the input VCF file

    :param infile: temporary input VCF file
    :param output: temporary output VCF file
    '''

    script = os.path.join(scripts_directory, 'flagPos.py')
    db = config['cosmic_params']['db']

    cmd = ['python', script, '--infile', infile, '--db', db, '--out', output,
           '--label', 'Cosmic', '--input_type', 'snv', '--flag_with_id', 'True']

    pypeliner.commandline.execute(*cmd)

#=======================================================================================================================
# Parse strelka vcf into csv format
#=======================================================================================================================

def parse_strelka(infile, output):
    parser = ParseStrelka(infile=infile, tid='NA', nid='NA', output=output,
                        keep_dbsnp=True,keep_1000gen=True,
                        remove_duplicates=True)

    parser.main()
