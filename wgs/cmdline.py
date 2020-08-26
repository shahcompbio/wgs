'''
Created on Feb 19, 2018

@author: dgrewal
'''
import json

import argparse
import pypeliner
from wgs import __version__


def add_global_args(subparser):
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
    alignment = add_global_args(alignment)

    # ================
    # realignment
    # ================
    realignment = subparsers.add_parser("realignment")
    realignment.set_defaults(which='realignment')
    realignment = add_global_args(realignment)

    # ================
    # variant calling
    # ================
    somatic_calling = subparsers.add_parser("somatic_calling")
    somatic_calling.set_defaults(which='somatic_calling')
    somatic_calling = add_global_args(somatic_calling)
    somatic_calling.add_argument(
        "--is_exome",
        default=False,
        action='store_true',
        help='''strelka filter'''
    )

    # ================
    # germline calling
    # ================
    germline_calling = subparsers.add_parser("germline_calling")
    germline_calling.set_defaults(which='germline_calling')
    add_global_args(germline_calling)

    # ================
    # breakpoints calling
    # ================
    sv_calling = subparsers.add_parser("breakpoint_calling")
    sv_calling.set_defaults(which='breakpoint_calling')
    sv_calling = add_global_args(sv_calling)

    # ================
    # copy number calling
    # ================
    cna_calling = subparsers.add_parser("copynumber_calling")
    cna_calling.set_defaults(which='copynumber_calling')
    cna_calling = add_global_args(cna_calling)
    cna_calling.add_argument(
        "--titan",
        default=False,
        action='store_true',
        help='''titan'''
    )
    cna_calling.add_argument(
        "--hmmcopy",
        default=False,
        action='store_true',
        help='''hmmcopy'''
    )
    cna_calling.add_argument(
        "--remixt",
        default=False,
        action='store_true',
        help='''remixt'''
    )

    # ================
    # postprocessing
    # ================
    postprocessing = subparsers.add_parser("postprocessing")
    postprocessing.set_defaults(which='postprocessing')
    postprocessing = add_global_args(postprocessing)
    postprocessing.add_argument(
        "--qc_metadata",
        required=True,
        help='''qc_metadata'''
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
