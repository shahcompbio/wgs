

def juno_config(reference, containers):
    docker_containers = containers['docker']

    if reference == 'grch37':
        reference = "/juno/work/shah/reference/genomes/gr37_decoys/gr37.fasta"
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
                'snpeff_config': '/juno/work/shah/reference/wgs_pipeline/snpEff.config'
            },
            'mutation_assessor_params': {
                'db': '/juno/work/shah/reference/databases/MA.hg19_v2/'
            },
            'dbsnp_params': {
                'db': '/juno/work/shah/reference/databases/dbsnp_142.human_9606.all.vcf.gz'
            },
            'thousandgen_params': {
                'db': '/juno/work/shah/reference/databases/1000G_release_20130502_genotypes.vcf.gz'
            },
            'cosmic_params': {
                'db': '/juno/work/shah/reference/databases/CosmicMutantExport.sorted.vcf.gz'
            },
            'mappability_ref': '/juno/work/shah/reference/wgs_pipeline/mask_regions_blacklist_crg_align36_table.txt',
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample': '/juno/work/shah/reference/wgs_pipeline/single_sample_plot_reference.h5',
            'thousandgen_params': {
                'db': '/juno/work/shah/reference/databases/1000G_release_20130502_genotypes.vcf.gz'
            },
            'dbsnp_params': {
                'db': '/juno/work/shah/reference/databases/dbsnp_142.human_9606.all.vcf.gz'
            },
        },
        'parse_strelka': {
            'keep_1000gen': True,
            ## TODO: why is this missing
            # 'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'filter_low_mappability': True,
        },
        'parse_museq': {
            'keep_1000gen': True,
            'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'pr_threshold': 0.85,
            'filter_low_mappability': True,
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
            'museqportrait': docker_containers['museqportrait'],
            'vizutils': docker_containers['vizutils'],
        }
    }

    sv_calling = {
        'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
        'samtools': 'samtools',
        'lumpyexpress': 'lumpyexpress',
        'refdata_destruct': '/juno/work/shah/reference/reference-grch37-decoys-destruct',
        'parse_lumpy': {
            'foldback_threshold': None,
            'mappability_ref': None,
            'chromosomes': map(str, range(1, 23) + ['X']),
            'normal_id': None,
            'deletion_size_threshold': 0,
            'tumour_read_support_threshold': 0,
            'project': None,
            'tumour_id': None,
            'confidence_interval_size': 500,
        },
        'parse_destruct': {
            'normal_id': None,
            'tumour_id': None,
            'case_id': None,
            'genes': None,
            'gene_locations': None,
            'chromosomes': map(str, range(1, 23) + ['X']),
            'deletion_size_threshold': 1000,
            'project': None,
            'types': None,
            'mappability_ref': '/juno/work/shah/reference/wgs_pipeline/mask_regions_blacklist_crg_align36_table_destruct.txt',
            'foldback_threshold': 30000,
            'readsupport_threshold': 4,
            'breakdistance_threshold': 30,
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
        'dbsnp_positions': '/juno/work/shah/reference/databases/common_all_dbSNP138.pos',
        'readcounter': {'w': 1000, 'q': 0},
        'correction': {
            'gc': '/juno/work/shah/reference/hmmcopy/GRCh37-lite.gc.ws_1000.wig',
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
        'pygenes_gtf': '/juno/work/shah/reference/wgs_pipeline/Homo_sapiens.GRCh37.73.gtf',
        'remixt_refdata': '/juno/work/shah/reference/reference-grch37-decoys-remixt',
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
        'titan_params': {
            'y_threshold': 20,
            'genome_type': 'NCBI',
            'map': '/juno/work/shah/reference/hmmcopy/GRCh37-lite.map.ws_1000.wig',
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
            'vizutils': docker_containers['vizutils'],
            'mutationseq': docker_containers['mutationseq'],
            'vcftools': docker_containers['vcftools'],
            'remixt': docker_containers['remixt']
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
            'SM': '{sample_id}',
            'PU': '{lane_id}',
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
