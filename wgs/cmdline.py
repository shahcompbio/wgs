'''
Created on Feb 19, 2018

@author: dgrewal
'''
import json

import argparse
import pypeliner

from wgs import __version__


def add_global_args(subparser, input_yaml=False):

    if input_yaml:
        subparser.add_argument("--input_yaml",
                               required=True,
                               help='''yaml file with tumour, normal and sampleids''')

    subparser.add_argument("--out_dir",
                           required=True,
                           help='''Path to output directory.''')

    subparser.add_argument("--config_file",
                           help='''Path to the config file.''')

    subparser.add_argument("--config_override",
                           type=json.loads,
                           help='''json string to override the defaults in config''')

    subparser.add_argument("--single_node",
                           default=False,
                           action='store_true',
                           help='''azure specific mode''')

    subparser.add_argument("--refdir",
                           required=True,
                           help='''reference data dir''')

    pypeliner.app.add_arguments(subparser)

    return subparser


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--version', action='version',
                        version='{version}'.format(version=__version__))

    subparsers = parser.add_subparsers()

    # ================
    # alignment
    # ================
    alignment = subparsers.add_parser("alignment")
    alignment.set_defaults(which='alignment')
    alignment = add_global_args(alignment, input_yaml=True)
    alignment.add_argument(
        "--picard_mem",
        default=8,
        type=int,
        help='''picard mem usage'''
    )

    # ================
    # realignment
    # ================
    realignment = subparsers.add_parser("realignment")
    realignment.set_defaults(which='realignment')
    realignment = add_global_args(realignment, input_yaml=True)
    realignment.add_argument(
        "--picard_mem",
        default=8,
        type=int,
        help='''picard mem usage'''
    )
    realignment.add_argument(
        "--ignore_bamtofastq_exception",
        default=False,
        action='store_true',
        help='''ignore the exception from bamtofastq'''
    )

    # ================
    # variant calling
    # ================
    somatic_calling = subparsers.add_parser("somatic_calling")
    somatic_calling.set_defaults(which='somatic_calling')
    somatic_calling = add_global_args(somatic_calling, input_yaml=True)
    somatic_calling.add_argument(
        "--is_exome",
        default=False,
        action='store_true',
        help='''strelka filter'''
    )

    # =========================
    # mutect - panel of normals
    # =========================
    somatic_panel_of_normals = subparsers.add_parser("somatic_panel_of_normals")
    somatic_panel_of_normals.set_defaults(which='somatic_panel_of_normals')
    somatic_panel_of_normals = add_global_args(somatic_panel_of_normals, input_yaml=True)

    # =========================
    # mutect
    # =========================
    mutect = subparsers.add_parser("mutect")
    mutect.set_defaults(which='mutect')
    mutect = add_global_args(mutect)
    mutect.add_argument(
        "--tumour",
        help='''strelka filter'''
    )
    mutect.add_argument(
        "--normal",
        help='''strelka filter'''
    )
    mutect.add_argument(
        "--tumour_id",
        help='''strelka filter'''
    )
    mutect.add_argument(
        "--normal_id",
        help='''strelka filter'''
    )

    # ================
    # germline calling
    # ================
    germline_calling = subparsers.add_parser("germline_calling")
    germline_calling.set_defaults(which='germline_calling')
    add_global_args(germline_calling, input_yaml=True)

    # ================
    # breakpoints calling
    # ================
    sv_calling = subparsers.add_parser("breakpoint_calling")
    sv_calling.set_defaults(which='breakpoint_calling')
    add_global_args(sv_calling, input_yaml=True)
    sv_calling.add_argument(
        "--svaba",
        default=False,
        action='store_true',
        help='''svaba'''
    )

    # ================
    # single sample copy number calling
    # ================
    single_sample_copynumber_calling = subparsers.add_parser("single_sample_copynumber_calling")
    single_sample_copynumber_calling.set_defaults(which='single_sample_copynumber_calling')
    add_global_args(single_sample_copynumber_calling, input_yaml=True)

    # ================
    # copy number calling
    # ================
    cna_calling = subparsers.add_parser("copynumber_calling")
    cna_calling.set_defaults(which='copynumber_calling')
    cna_calling = add_global_args(cna_calling, input_yaml=True)
    cna_calling.add_argument(
        "--titan",
        default=False,
        action='store_true',
        help='''titan'''
    )
    cna_calling.add_argument(
        "--remixt",
        default=False,
        action='store_true',
        help='''remixt'''
    )

    # ================
    # sample_qc
    # ================
    sample_qc = subparsers.add_parser("sample_qc")
    sample_qc.set_defaults(which='sample_qc')
    sample_qc = add_global_args(sample_qc, input_yaml=True)
    sample_qc.add_argument(
        '--bins',
        default=2000
    )
    sample_qc.add_argument(
        '--mapping_qual_threshold',
        default=0
    )
    sample_qc.add_argument(
        '--normal_only',
        action='store_true',
        default=False
    )

    # ================
    # cohort qc
    # ================
    cohort_qc = subparsers.add_parser("cohort_qc")
    cohort_qc.set_defaults(which='cohort_qc')
    cohort_qc = add_global_args(cohort_qc, input_yaml=True)
    cohort_qc.add_argument(
        "--API_key",
        default=False,
        help='''API key for account on oncokv to run MafAnnotate.py. '''
    )

    # ======================================
    # generates pipeline and batch configs
    # ======================================
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("--pipeline_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--batch_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--config_override",
                                 type=json.loads,
                                 help='''json string to override the defaults in config''')

    args = vars(parser.parse_args())

    return args
