"""
Created on Feb 19, 2018

@author: dgrewal
"""
from alignment import alignment_workflow
from cmdline import parse_args
from copynumber_calling import copynumber_calling_workflow
from realign import realign_bam_workflow
from breakpoint_calling import breakpoint_calling_workflow
from variant_calling import variant_calling_workflow
from postprocessing import postprocessing_workflow
from wgs.config import batch_config


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

    if args["which"] == "alignment":
        args = generate_config(args)
        alignment_workflow(args)

    if args["which"] == "variant_calling":
        args = generate_config(args)
        variant_calling_workflow(args)

    if args["which"] == "breakpoint_calling":
        args = generate_config(args)
        breakpoint_calling_workflow(args)

    if args["which"] == "copynumber_calling":
        args = generate_config(args)
        copynumber_calling_workflow(args)

    if args["which"] == "realignment":
        args = generate_config(args)
        realign_bam_workflow(args)

    if args["which"] == "postprocessing":
        args = generate_config(args)
        postprocessing_workflow(args)



if __name__ == "__main__":
    main()
