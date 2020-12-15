'''
Created on Feb 21, 2018

@author: dgrewal
'''
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.config import config


def create_somatic_consensus_workflow(
        mutect_snv_vcf,
        strelka_snv_vcf,
        strelka_indel_vcf,
        museq_snv_vcf,
        consensus_maf,
        chromosomes,
        reference_vep,
        normal_id,
        tumour_id,
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='snv_consensus',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.somatic_calling_consensus.consensus.main',
        args=(
            mgd.InputFile(museq_snv_vcf),
            mgd.InputFile(strelka_snv_vcf),
            mgd.InputFile(mutect_snv_vcf),
            mgd.InputFile(strelka_indel_vcf),
            mgd.TempOutputFile('consensus.vcf'),
            mgd.TempOutputFile('counts.csv'),
            chromosomes,
        ),
    )

    workflow.subworkflow(
        name="consensus_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.TempInputFile('consensus.vcf'),
            mgd.TempOutputFile('consensus.maf'),
            reference_vep,
        ),
        kwargs={'normal_id': normal_id, 'tumour_id': tumour_id}
    )

    workflow.transform(
        name='maf_counts',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.somatic_calling_consensus.tasks.update_maf_counts',
        args=(
            mgd.TempInputFile('consensus.maf'),
            mgd.TempInputFile('counts.csv'),
            mgd.OutputFile(consensus_maf),
        )
    )

    return workflow
