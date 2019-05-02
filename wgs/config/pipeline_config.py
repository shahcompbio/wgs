import os
import warnings
import yaml
from wgs.utils import helpers


def luna_config(reference):

    if reference == 'grch37':
        reference = "/juno/work/shah/reference/gr37.fasta"
    else:
        reference = None

    globals = {
        'memory': {'low': 5, 'med': 10, 'high': 15,},
        'threads': 1,
    }

    variant_calling = {
        'split_size': 1e7,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'reference': reference,
        'annotation_params': {
            'snpeff_params': {
                'snpeff_config': '/juno/work/shah/reference/snpEff.config'
            },
            'mutation_assessor_params': {
                'db': '/juno/work/shah/reference/MA.hg19_v2/'
            },
            'dbsnp_params': {
                'db': '/juno/work/shah/reference/dbsnp_142.human_9606.all.vcf.gz'
             },
            'thousandgen_params': {
                'db': '/juno/work/shah/reference/1000G_release_20130502_genotypes.vcf.gz'
            },
            'cosmic_params': {
                'db': '/juno/work/shah/reference/CosmicMutantExport.sorted.vcf.gz'
            },
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample':'/juno/work/shah/reference/single_sample_plot_reference.h5',
            'thousandgen_params': {
                'db': '/juno/work/shah/reference/1000G_release_20130502_genotypes.vcf.gz'
            },
            'dbsnp_params': {
                'db': '/juno/work/shah/reference/dbsnp_142.human_9606.all.vcf.gz'
            },
        },
        'parse_strelka':{
            'keep_1000gen': True,
            ## TODO: why is this missing
            # 'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/juno/work/shah/reference/mask_regions_blacklist_crg_align36_table.txt',
         },
        'parse_museq': {
            'keep_1000gen': True,
            'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/juno/work/shah/reference/mask_regions_blacklist_crg_align36_table.txt',
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
        }
    }

    sv_calling = {
        'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
        'samtools': 'samtools',
        'lumpyexpress': 'lumpyexpress',
        'refdata_destruct': '/juno/work/shah/reference/',
        'parse_lumpy':{
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
            'mappability_ref': '/juno/work/shah/reference/mask_regions_blacklist_crg_align36_table_destruct.txt',
            'foldback_threshold': 30000,
            'readsupport_threshold': 4,
            'breakdistance_threshold': 30,
        },
    }

    cna_calling = {
        'split_size': 1e7,
        "min_num_reads": 5,
        "reference_genome": reference,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'dbsnp_positions': '/juno/work/shah/reference/common_all_dbSNP138.pos',
        'readcounter': {'w': 1000, 'q': 0},
        'correction': {
            'gc': '/juno/work/shah/reference/GRCh37-lite.gc.ws_1000.wig',
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
        'pygenes_gtf': '/juno/work/shah/reference/Homo_sapiens.GRCh37.73.gtf',
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
            'map': '/juno/work/shah/reference/GRCh37-lite.map.ws_1000.wig',
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
        }
    }

    alignment = {
        "ref_genome":{
            'file': reference,
            'header':{
                'UR': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome',
                'AS': 'hg19/1000genomes',
                'SP': 'Homo sapiens'
            }
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
        }
    }

    config = locals()

    return config


def azure_config(reference):

    if reference == 'grch37':
        reference = "/refdata/GRCh37-lite.fa"
    else:
        reference = None

    globals = {
        'memory': {'low': 5, 'med': 10, 'high': 15,},
        'threads': 1,
    }

    variant_calling = {
        'split_size': 1e7,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'reference': reference,
        'annotation_params': {
            'snpeff_params': {
                'snpeff_config': '/refdata/snpEff.config'
            },
            'mutation_assessor_params': {
                'db': '/refdata/MA.hg19_v2/'
            },
            'dbsnp_params': {
                'db': '/refdata/dbsnp_142.human_9606.all.vcf.gz'
            },
            'thousandgen_params': {
                'db': '/refdata/1000G_release_20130502_genotypes.vcf.gz'
            },
            'cosmic_params': {
                'db': '/refdata/CosmicMutantExport.sorted.vcf.gz'
            },
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample': '/refdata/single_sample_plot_reference.h5',
            'thousandgen_params': {
                'db': '/refdata/1000G_release_20130502_genotypes.vcf.gz'
            },
            'dbsnp_params': {
                'db': '/refdata/dbsnp_142.human_9606.all.vcf.gz'
            },
        },
        'parse_strelka':{
            'keep_1000gen': True,
            ## TODO: why is this missing
            # 'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/refdata/mask_regions_blacklist_crg_align36_table.txt',
         },
        'parse_museq': {
            'keep_1000gen': True,
            'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/refdata/mask_regions_blacklist_crg_align36_table.txt',
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
        }
    }

    sv_calling = {
        'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
        'samtools': 'samtools',
        'lumpyexpress': 'lumpyexpress',
        'refdata_destruct': '/refdata/',
        'parse_lumpy':{
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
            'mappability_ref': '/refdata/mask_regions_blacklist_crg_align36_table_destruct.txt',
            'foldback_threshold': 30000,
            'readsupport_threshold': 4,
            'breakdistance_threshold': 30
        }
    }

    cna_calling = {
        'split_size': 1e7,
        "min_num_reads": 5,
        "reference_genome": reference,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'dbsnp_positions': '/refdata/common_all_dbSNP138.pos',
        'readcounter': {'w': 1000, 'q': 0},
        'correction' : {
            'gc': '/refdata/GRCh37-lite.gc.ws_1000.wig',
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
        'pygenes_gtf': '/refdata/Homo_sapiens.GRCh37.73.gtf',
        'remixt_refdata': '/refdata/reference-grch37-decoys-remixt',
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
            'map': '/refdata/GRCh37-lite.map.ws_1000.wig',
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
        }
    }

    alignment = {
        "ref_genome":{
            'file': reference,
            'header':{
                'UR': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome',
                'AS': 'hg19/1000genomes',
                'SP': 'Homo sapiens'
            }
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
        }
    }


    config = locals()

    return config



def shahlab_config(reference):

    if reference == 'grch37':
        reference = "/shahlab/pipelines/reference/GRCh37-lite.fa"
    else:
        reference = None

    globals = {
        'memory': {'low': 5, 'med': 10, 'high': 15,},
        'threads': 1,
    }

    variant_calling = {
        'split_size': 1e7,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'reference': reference,
        'annotation_params': {
            'snpeff_params': {
                'snpeff_config': '/shahlab/pipelines/reference/snpEff.config'
            },
            'mutation_assessor_params': {
                'db': '/shahlab/pipelines/reference/MA.hg19_v2/'
            },
            'dbsnp_params': {
                'db': '/shahlab/pipelines/reference/dbsnp_142.human_9606.all.vcf.gz'
            },
            'thousandgen_params': {
                'db': '/shahlab/pipelines/reference/1000G_release_20130502_genotypes.vcf.gz'
            },
            'cosmic_params': {
                'db': '/shahlab/pipelines/reference/CosmicMutantExport.sorted.vcf.gz'
            },
        },
        'plot_params': {
            'threshold': 0.5,
            'refdata_single_sample': '/shahlab/pipelines/reference/single_sample_plot_reference.h5',
            'thousandgen_params': {
                'db': '/shahlab/pipelines/reference/1000G_release_20130502_genotypes.vcf.gz'
            },
            'dbsnp_params': {
                'db': '/shahlab/pipelines/reference/dbsnp_142.human_9606.all.vcf.gz'
            },
        },
        'parse_strelka': {
            'keep_1000gen': True,
            ## TODO: why is this missing
            # 'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/shahlab/pipelines/reference/mask_regions_blacklist_crg_align36_table.txt',
        },
        'parse_museq': {
            'keep_1000gen': True,
            'keep_cosmic': True,
            'remove_duplicates': False,
            'keep_dbsnp': True,
            'chromosomes': map(str, range(1, 23)) + ['X'],
            'mappability_ref': '/shahlab/pipelines/reference/mask_regions_blacklist_crg_align36_table.txt',
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
        }
    }

    sv_calling = {
        'extractSplitReads_BwaMem': 'lumpy_extractSplitReads_BwaMem',
        'samtools': 'samtools',
        'lumpyexpress': 'lumpyexpress',
        'refdata_destruct': '/shahlab/pipelines/reference/refdir_destruct_GRCh37/',
        'parse_lumpy':{
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
            'mappability_ref': '/shahlab/pipelines/reference/mask_regions_blacklist_crg_align36_table_destruct.txt',
            'foldback_threshold': 30000,
            'readsupport_threshold': 4,
            'breakdistance_threshold': 30
        }
    }

    cna_calling = {
        'split_size': 1e7,
        "min_num_reads": 5,
        "reference_genome": reference,
        'chromosomes': map(str, range(1, 23) + ['X']),
        'dbsnp_positions': '/shahlab/pipelines/reference/common_all_dbSNP138.pos',
        'readcounter': {'w': 1000, 'q': 0},
        'correction': {
            'gc': '/shahlab/pipelines/reference/GRCh37-lite.gc.ws_1000.wig',
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
        'pygenes_gtf': '/shahlab/pipelines/reference/Homo_sapiens.GRCh37.73.gtf',
        'remixt_refdata': '/shahlab/pipelines/reference/refdir_remixt_new_GRCh37/',
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
            'map': '/shahlab/pipelines/reference/GRCh37-lite.map.ws_1000.wig',
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
        }
    }

    alignment = {
        "ref_genome":{
            'file': reference,
            'header':{
                'UR': 'http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome',
                'AS': 'hg19/1000genomes',
                'SP': 'Homo sapiens'
            }
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
        }
    }

    config = locals()

    return config

def get_config(override):

    if override["cluster"] == "shahlab":
        config = shahlab_config(override["reference"])
    elif override["cluster"] == "juno":
        config = luna_config(override["reference"])
    elif override["cluster"] == "azure":
        config = azure_config(override["reference"])
    else:
        raise Exception()

    return config


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)


def generate_pipeline_config(args):

    if args['which'] == 'generate_config':
        config_yaml = args['pipeline_config']
        config_yaml = os.path.abspath(config_yaml)
    else:
        config_yaml = "config.yaml"
        tmpdir = args.get("tmpdir", None)
        pipelinedir = args.get("pipelinedir", None)

        # use pypeliner tmpdir to store yaml
        if pipelinedir:
            config_yaml = os.path.join(pipelinedir, config_yaml)
        elif tmpdir:
            config_yaml = os.path.join(tmpdir, config_yaml)
        else:
            warnings.warn("no tmpdir specified, generating configs in working dir")
            config_yaml = os.path.join(os.getcwd(), config_yaml)

        config_yaml = helpers.get_incrementing_filename(config_yaml)
    print config_yaml

    params_override = {'cluster': 'azure', 'reference': 'grch37'}
    if args['config_override']:
        params_override.update(args["config_override"])

    helpers.makedirs(config_yaml, isfile=True)

    config = get_config(params_override)
    write_config(config, config_yaml)

    args["config_file"] = config_yaml

    print config_yaml
    return args

