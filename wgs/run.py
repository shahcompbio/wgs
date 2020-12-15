"""
Created on Feb 19, 2018

@author: dgrewal
"""
from wgs.alignment import alignment_workflow
from wgs.breakpoint_calling import breakpoint_calling_workflow
from wgs.cmdline import parse_args
from wgs.cohort_qc import cohort_qc_workflow
from wgs.config import batch_config
from wgs.copynumber_calling import copynumber_calling_workflow
from wgs.germline_calling import germline_calling_workflow
from wgs.realign import realign_bam_workflow
from wgs.sample_qc import sample_qc_workflow
from wgs.somatic_calling import somatic_calling_workflow


def generate_config(args):
    if args['which'] == 'generate_config':
        if args.get("batch_config"):
            args = batch_config.generate_submit_config_in_temp(args)
    else:
        if not args.get("submit_config"):
            args = batch_config.generate_submit_config_in_temp(args)
    return args


def main():
    args = parse_args()

    if args["which"] == "generate_config":
        generate_config(args)

    if args["which"] == "cohort_qc":
        args = generate_config(args)
        cohort_qc_workflow(args)

    if args["which"] == "alignment":
        args = generate_config(args)
        alignment_workflow(args)

    if args["which"] == "somatic_calling":
        args = generate_config(args)
        somatic_calling_workflow(args)

    if args["which"] == "germline_calling":
        args = generate_config(args)
        germline_calling_workflow(args)

    if args["which"] == "breakpoint_calling":
        args = generate_config(args)
        breakpoint_calling_workflow(args)

    if args["which"] == "copynumber_calling":
        args = generate_config(args)
        copynumber_calling_workflow(args)

    if args["which"] == "realignment":
        args = generate_config(args)
        realign_bam_workflow(args)

    if args["which"] == "sample_qc":
        args = generate_config(args)
        sample_qc_workflow(args)


if __name__ == "__main__":
    main()
