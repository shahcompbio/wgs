"""
Created on Feb 19, 2018

@author: dgrewal
"""
import os
import pypeliner
import pypeliner.managed as mgd
from cmdline import parse_args
from variant_calling import variant_calling_workflow
from sv_calling import sv_calling_workflow
from cna_calling import cna_calling_workflow
from wgs.config import pipeline_config
from wgs.config import batch_config

def generate_config(args):
    if args['which'] == 'generate_config':
        if args.get("pipeline_config"):
            args = pipeline_config.generate_pipeline_config(args)
        if args.get("batch_config"):
            args = batch_config.generate_pipeline_config(args)
    else:
        if not args.get("config_file"):
            args = pipeline_config.generate_pipeline_config(args)
        if not args.get("submit_config"):
            args = batch_config.generate_batch_config(args)
    return args


def main():
    args = parse_args()

    if args["which"] == "generate_config":
        generate_config(args)

    if args["which"] == "variant_calling":
        args = generate_config(args)
        variant_calling_workflow(args)

    if args["which"] == "sv_calling":
        args = generate_config(args)
        sv_calling_workflow(args)

    if args["which"] == "cna_calling":
        args = generate_config(args)
        cna_calling_workflow(args)




if __name__ == "__main__":
    main()


