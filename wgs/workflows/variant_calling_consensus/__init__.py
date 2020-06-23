'''
Created on Feb 21, 2018

@author: dgrewal
'''
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.config import config


def create_consensus_workflow(
        museq_snv,
        strelka_snv,
        strelka_indel,
        somatic_calls,
        somatic_snpeff,
        somatic_ma,
        somatic_ids,
        indel_calls,
        indel_snpeff,
        indel_ma,
        indel_ids,
        refdir
):
    params = config.default_params('variant_calling')
    chromosomes = config.refdir_data(refdir)['params']['chromosomes']

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_strelka_indel',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_indel, extensions=['.csi', '.tbi']),
            mgd.OutputFile(indel_calls, extensions=['.yaml']),
            mgd.OutputFile(indel_snpeff, extensions=['.yaml']),
            mgd.OutputFile(indel_ma, extensions=['.yaml']),
            mgd.OutputFile(indel_ids, extensions=['.yaml']),
            params["parse_strelka"],
            chromosomes,
            mgd.TempSpace("tempdir_strelka_indel")
        ),
    )

    workflow.transform(
        name='parse_museq_snv',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(museq_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('museq_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ids.csv', extensions=['.yaml']),
            params["parse_museq"],
            chromosomes,
            mgd.TempSpace("tempdir_parse_museq_snv")
        ),
    )

    workflow.transform(
        name='parse_strelka_snv',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('strelka_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ids.csv', extensions=['.yaml']),
            params["parse_strelka"],
            chromosomes,
            mgd.TempSpace("tempdir_parse_strelka_snv")
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_snv.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_calls, extensions=['.yaml']),
        ),
    )

    workflow.transform(
        name='merge_snpeff',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_snpeff.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_snpeff, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ma',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ma.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ma.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ma, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ids',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ids.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ids.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ids, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    return workflow
