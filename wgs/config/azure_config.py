def azure_config(reference, containers):
    docker_containers = containers['docker']

    if reference == 'grch37':
        reference = " /refdata/GRCh37-lite/GRCh37-lite.fa"
    else:
        reference = None

    globals = {
        'memory': {'low': 5, 'med': 10, 'high': 15, },
        'threads': 8,
    }

    variant_calling = {
        'split_size': 1e7,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'reference': reference,
        'annotation_params': {
            'snpeff_params': {
                'snpeff_config': '/refdata/wgs_pipeline/snpEff.config'
            },
            'mutation_assessor_params': {
                'db': '/refdata/databases/MA.hg19_v2/'
            },
            'dbsnp_params': {
                'db': '/refdata/databases/dbsnp_142.human_9606.all.vcf.gz'
            },
            'thousandgen_params': {
                'db': '/refdata/databases/1000G_release_20130502_genotypes.vcf.gz'
            },
            'cosmic_params': {
                'db': '/refdata/databases/CosmicMutantExport.sorted.vcf.gz'
            },
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample': '/refdata/wgs_pipeline/single_sample_plot_reference.h5',
            'thousandgen_params': {
                'db': '/refdata/databases/1000G_release_20130502_genotypes.vcf.gz'
            },
            'dbsnp_params': {
                'db': '/refdata/databases/dbsnp_142.human_9606.all.vcf.gz'
            },
        },
        'parse_strelka': {
            'keep_1000gen': True,
            ## TODO: why is this missing
            # 'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt',
        },
        'parse_museq': {
            'keep_1000gen': True,
            'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt',
            'pr_threshold': 0.85
        },
        'museq_params': {
            'threshold': 0.5,
            'verbose': True,
            'purity': 70,
            'coverage': 4,
            'buffer_size': '2G',
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
            'mutationseq': docker_containers['mutationseq'],
            'vizutils': docker_containers['vizutils'],
        }
    }

    sv_calling = {
        'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
        'samtools': 'samtools',
        'lumpyexpress': 'lumpyexpress',
        'refdata_destruct': '/refdata/human/',
        'destruct_config': {
            'genome_fasta': '/refdata/human/GRCh37-lite.fa',
            'genome_fai': '/refdata/human/GRCh37-lite.fa.fai',
            'gtf_filename': '/refdata/human/GRCh37-lite.gtf',
        },
        'mappability_ref': '/refdata/wgs_pipeline/mask_regions_blacklist_crg_align36_table_destruct.txt',
        'parse_lumpy': {
            'chromosomes': map(str, range(1, 23) + ['X']),
            'deletion_size_threshold': 0,
            'tumour_read_support_threshold': 0,
        },
        'parse_destruct': {
            'chromosomes': map(str, range(1, 23) + ['X']),
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
            'vizutils': docker_containers['vizutils'],
        }
    }

    copynumber_calling = {
        'ncpus': 8,
        'split_size': 1e7,
        "min_num_reads": 5,
        "reference_genome": reference,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'dbsnp_positions': '/refdata/databases/common_all_dbSNP138.pos',
        'readcounter': {'w': 1000, 'q': 0},
        'correction': {
            'gc': '/refdata/hmmcopy/GRCh37-lite.gc.ws_1000.wig',
            'map': '/refdata/hmmcopy/GRCh37-lite.map.ws_1000.wig',
            'map_cutoff': 0.85
        },
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
        'pygenes_gtf': '/refdata/wgs_pipeline/Homo_sapiens.GRCh37.73.gtf',
        'remixt_refdata': '/refdata/reference-remixt',
        'museq_params': {
            'threshold': 0.85,
            'verbose': True,
            'purity': 70,
            'coverage': 4,
            'buffer_size': '2G',
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
            'map': '/refdata/hmmcopy/GRCh37-lite.map.ws_1000.wig',
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
            'chromosomes': map(str, range(1, 23) + ['X']),
            'genes': None,
            'types': None,
        },
        'docker': {
            'wgs': docker_containers['wgs'],
            'titan': docker_containers['titan'],
            'hmmcopy': docker_containers['hmmcopy'],
            'vizutils': docker_containers['vizutils'],
            'mutationseq': docker_containers['mutationseq'],
            'vcftools': docker_containers['vcftools'],
        }
    }

    alignment = {
        "ref_genome": {
            'file': reference,
            'header': {
                'UR': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome',
                'AS': 'hg19/1000genomes',
                'SP': 'Homo sapiens'
            }
        },
        'picard_wgs_params': {
            "min_bqual": 20,
            "min_mqual": 20,
            "count_unpaired": False,
        },
        'threads': 8,
        'aligner': 'bwa-mem',
        'split_size': 1e7,
        'read_group_info': {
            'ID': '{sample_id}_{lane_id}',
            'PU': '{lane_id}',
            'SM': '{sample_id}',
            'CN': 'IGO_MSKCC',
            'PL': 'ILLUMINA',
        },
        'docker': {
            'wgs': docker_containers['wgs'],
            'bwa': docker_containers['bwa'],
            'samtools': docker_containers['samtools'],
            'picard': docker_containers['picard'],
            'fastqc': docker_containers['fastqc']
        }
    }

    config = locals()

    return config
