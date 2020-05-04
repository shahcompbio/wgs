import os

import wgs
import yaml


def default_params(mode='all'):
    variant_calling = {
        'split_size': 1e7,
        'strelka_depth_threshold': True,
        'germline_portrait_threshold': 0.5,
        'parse_strelka': {},
        'parse_museq': {
            'pr_threshold': 0.85
        },
        'museq_params': {
            'threshold': 0.5,
            'verbose': True,
            'coverage': 4,
            'mapq_threshold': 10,
            'indl_threshold': 0.05,
            'normal_variant': 25,
            'tumour_variant': 2,
            'baseq_threshold': 20,
        },
    }

    breakpoint_calling = {
        'chromosomes': list(map(str, range(1, 23))) + ['X'],
        'lumpy_paths': {
            'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
            'samtools': 'samtools',
            'lumpyexpress': 'lumpyexpress',
        },
        'parse_lumpy': {
            'deletion_size_threshold': 0,
            'tumour_read_support_threshold': 0,
        },
        'parse_destruct': {
            'deletion_size_threshold': 1000,
            'foldback_threshold': 30000,
            'readsupport_threshold': 4,
            'breakdistance_threshold': 30
        },
        'consensus': {
            'confidence_interval_size': 500,
        },
    }

    copynumber_calling = {
        'split_size': 1e7,
        'chromosomes': list(map(str, range(1, 23))) + ['X'],
        'readcounter': {'w': 1000, 'q': 0},
        'map_cutoff': 0.85,
        'genome_type': 'NCBI',
        'titan_intervals': [
            {'num_clusters': 1, 'ploidy': 2},
            {'num_clusters': 2, 'ploidy': 2},
            {'num_clusters': 3, 'ploidy': 2},
            {'num_clusters': 4, 'ploidy': 2},
            {'num_clusters': 5, 'ploidy': 2},
            {'num_clusters': 1, 'ploidy': 4},
            {'num_clusters': 2, 'ploidy': 4},
            {'num_clusters': 3, 'ploidy': 4},
            {'num_clusters': 4, 'ploidy': 4},
            {'num_clusters': 5, 'ploidy': 4},
        ],
        'museq_params': {
            'threshold': 0.85,
            'verbose': True,
            'mapq_threshold': 10,
            'indl_threshold': 0.05,
            'normal_variant': 25,
            'tumour_variant': 2,
            'baseq_threshold': 10,
        },
        'hmmcopy_params': {
            'normal_copy': None,
            'normal_table': None,
            'normal_table_out': None,
            'm': None,
            'mu': None,
            'kappa': None,
            'e': None,
            'S': None,
            'strength': None,
            'lambda': None,
            'nu': None,
            'eta': None,
            'gamma': None,
        },
        'titan_params': {
            'y_threshold': 20,
            'num_cores': 4,
            'myskew': 0,
            'estimate_ploidy': 'TRUE',
            'normal_param_nzero': 0.5,
            'normal_estimate_method': 'map',
            'max_iters': 50,
            'pseudo_counts': 1e-300,
            'txn_exp_len': 1e16,
            'txn_z_strength': 1e6,
            'alpha_k': 15000,
            'alpha_high': 20000,
            'max_copynumber': 8,
            'symmetric': 'TRUE',
            'chrom': 'NULL',
            'max_depth': 1000,
        },
        'parse_titan': {
            'segment_size_threshold': 5000,
            'genes': None,
            'types': None,
        },
    }

    alignment = {
        'picard_wgs_params': {
            "min_bqual": 20,
            "min_mqual": 20,
            "count_unpaired": False,
        },
        'threads': 8,
        'aligner': 'bwa-mem',
        'split_size': 1e7,
    }

    config = {
        'copynumber_calling': copynumber_calling,
        'breakpoint_calling': breakpoint_calling,
        'variant_calling': variant_calling,
        'alignment': alignment
    }

    if not mode=='all':
        return config[mode]

    return config


def refdir_data(refdir):
    yamlpath = os.path.join(refdir, 'metadata.yaml')
    with open(yamlpath) as yamlfile:
        yamldata = yaml.safe_load(yamlfile)

    for k, v in yamldata['paths'].items():
        yamldata['paths'][k] = os.path.join(refdir, v)

    return yamldata


def containers(container_name):
    version = wgs.__version__
    # strip setuptools metadata
    version = version.split("+")[0]

    docker_images = {
        'bwa': 'bwa:v0.0.1',
        'samtools': 'samtools:v0.0.2',
        'picard': 'picard:v0.0.1',
        'wgs': 'wgs:v{}'.format(version),
        'strelka': 'strelka:v0.0.1',
        'mutationseq': 'mutationseq:v0.0.1',
        'vcftools': 'vcftools:v0.0.1',
        'snpeff': 'vcftools:v0.0.1',
        'titan': 'titan:v0.0.2',
        'destruct': 'destruct:v0.0.2',
        'lumpy': 'lumpy:v0.0.1',
        'fastqc': 'fastqc:v0.0.1',
        'hmmcopy': 'hmmcopy:v0.0.1',
        'circos': 'circos:v0.0.1',
        'igvtools': 'igvtools:v0.0.1',
    }

    return docker_images[container_name]
