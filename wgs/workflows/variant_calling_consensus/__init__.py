'''
Created on Feb 21, 2018

@author: pwalters
'''
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_consensus_workflow(
        museq_germline,
        museq_snv,
        strelka_snv,
        strelka_indel,
        somatic_calls,
        indel_calls,
        germline_calls,
        outdir,
        global_config,
        varcall_config
):
    germline_snpeff_annotations = os.path.join(outdir, 'germline_snpeff_annotations.csv.gz')
    indel_snpeff_annotations = os.path.join(outdir, 'indel_snpeff_annotations.csv.gz')
    somatic_snpeff_annotations = os.path.join(outdir, 'somatic_snpeff_annotations.csv.gz')

    germline_ma_annotations = os.path.join(outdir, 'germline_ma_annotations.csv.gz')
    indel_ma_annotations = os.path.join(outdir, 'indel_ma_annotations.csv.gz')
    somatic_ma_annotations = os.path.join(outdir, 'somatic_ma_annotations.csv.gz')

    germline_ids_annotations = os.path.join(outdir, 'germline_ids_annotations.csv.gz')
    indel_ids_annotations = os.path.join(outdir, 'indel_ids_annotations.csv.gz')
    somatic_ids_annotations = os.path.join(outdir, 'somatic_ids_annotations.csv.gz')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='parse_museq_germlines',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(museq_germline, extensions=['.csi', '.tbi']),
            mgd.OutputFile(germline_calls, extensions=['.yaml']),
            mgd.OutputFile(germline_snpeff_annotations, extensions=['.yaml']),
            mgd.OutputFile(germline_ma_annotations, extensions=['.yaml']),
            mgd.OutputFile(germline_ids_annotations, extensions=['.yaml']),
            varcall_config["parse_museq"],
            mgd.TempSpace("tempdir_parse_germlines")
        ),
    )

    workflow.transform(
        name='parse_strelka_indel',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_indel, extensions=['.csi', '.tbi']),
            mgd.OutputFile(indel_calls, extensions=['.yaml']),
            mgd.OutputFile(indel_snpeff_annotations, extensions=['.yaml']),
            mgd.OutputFile(indel_ma_annotations, extensions=['.yaml']),
            mgd.OutputFile(indel_ids_annotations, extensions=['.yaml']),
            varcall_config["parse_strelka"],
            mgd.TempSpace("tempdir_strelka_indel")
        ),
    )

    workflow.transform(
        name='parse_museq_snv',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(museq_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('museq_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('museq_ids.csv', extensions=['.yaml']),
            varcall_config["parse_museq"],
            mgd.TempSpace("tempdir_parse_museq_snv")
        ),
    )

    workflow.transform(
        name='parse_strelka_snv',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.parse_vcf',
        args=(
            mgd.InputFile(strelka_snv, extensions=['.csi', '.tbi']),
            mgd.TempOutputFile('strelka_snv.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ma.csv', extensions=['.yaml']),
            mgd.TempOutputFile('strelka_snv_ids.csv', extensions=['.yaml']),
            varcall_config["parse_strelka"],
            mgd.TempSpace("tempdir_parse_strelka_snv")
        ),
    )

    workflow.transform(
        name='merge_snvs',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
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
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_snpeff.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_snpeff.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_snpeff_annotations, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ma',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ma.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ma.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ma_annotations, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    workflow.transform(
        name='merge_ids',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.variant_calling_consensus.tasks.merge_overlap',
        args=(
            [mgd.TempInputFile('strelka_snv_ids.csv', extensions=['.yaml']),
             mgd.TempInputFile('museq_ids.csv', extensions=['.yaml'])],
            mgd.OutputFile(somatic_ids_annotations, extensions=['.yaml']),
        ),
        kwargs={'on': ['chrom', 'pos']}
    )

    return workflow
