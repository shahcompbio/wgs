import os

import yaml


def load_refdir_metadata(refdir):
    yamlpath = os.path.join(refdir, 'metadata.yaml')

    with open(yamlpath) as yamlfile:
        yamldata = yaml.safe_load(yamlfile)

    return yamldata


def pipeline_config(containers, refdir):
    def get_path(metadata_key):
        return os.path.join(refdir, metadata[metadata_key])

    metadata = load_refdir_metadata(refdir)

    docker_containers = containers['docker']

    globals = {
        'memory': {'low': 5, 'med': 10, 'high': 15, },
        'threads': 8,
    }

    variant_calling = {
        'split_size': 1e7,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'reference': get_path('reference'),
        'strelka_depth_threshold': True,
        'databases': {
            'snpeff_params': {
                'snpeff_config': get_path('snpeff_config')
            },
            'mutation_assessor_params': {
                'db': get_path('mutation_assessor')
            },
            'dbsnp_params': {
                'db': get_path('dbsnp'),
            },
            'thousandgen_params': {
                'db': get_path('thousand_genomes')
            },
            'cosmic_params': {
                'db': get_path('cosmic')
            },
            'mappability_ref': get_path('blacklist')
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample': get_path('germline_portrait_ref'),
        },
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
        'docker': {
            'wgs': docker_containers['wgs'],
            'strelka': docker_containers['strelka'],
            'vcftools': docker_containers['vcftools'],
            'samtools': docker_containers['samtools'],
            'mutationseq': docker_containers['mutationseq'],
        }
    }

    sv_calling = {
        'chromosomes': map(str, range(1, 23) + ['X']),
        'lumpy_paths': {
            'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
            'samtools': 'samtools',
            'lumpyexpress': 'lumpyexpress',
        },
        'refdata_destruct': get_path('refdata_destruct'),
        'destruct_config': {
            'genome_fasta': get_path('reference'),
            'genome_fai': get_path('reference') + '.fai',
            'gtf_filename': get_path('gtf'),
        },
        'mappability_ref': get_path('blacklist_destruct'),
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
        'docker': {
            'wgs': docker_containers['wgs'],
            'destruct': docker_containers['destruct'],
            'lumpy': docker_containers['lumpy'],
            'samtools': docker_containers['samtools'],
        }
    }

    copynumber_calling = {
        'split_size': 1e7,
        "reference_genome": get_path('reference'),
        'chromosomes': map(str, range(1, 23) + ['X']),
        'dbsnp_positions': get_path('het_positions_titan'),
        'readcounter': {'w': 1000, 'q': 0},
        'reference_wigs': {
            'gc': get_path('gc_wig'),
            'map': get_path('map_wig'),
        },
        'map_cutoff': 0.85,
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
        'pygenes_gtf': get_path('gtf'),
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
            'genome_type': 'NCBI',
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
        'docker': {
            'wgs': docker_containers['wgs'],
            'titan': docker_containers['titan'],
            'hmmcopy': docker_containers['hmmcopy'],
            'mutationseq': docker_containers['mutationseq'],
            'vcftools': docker_containers['vcftools'],
        }
    }

    alignment = {
        "ref_genome": get_path('reference'),
        'picard_wgs_params': {
            "min_bqual": 20,
            "min_mqual": 20,
            "count_unpaired": False,
        },
        'threads': 8,
        'aligner': 'bwa-mem',
        'split_size': 1e7,
        'docker': {
            'wgs': docker_containers['wgs'],
            'bwa': docker_containers['bwa'],
            'samtools': docker_containers['samtools'],
            'picard': docker_containers['picard'],
            'fastqc': docker_containers['fastqc']
        }
    }

    config = {
        'copynumber_calling': copynumber_calling,
        'globals': globals,
        'sv_calling': sv_calling,
        'variant_calling': variant_calling,
        'alignment': alignment
    }

    return config
