'''
Created on Feb 21, 2018

@author: dgrewal
'''

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_germline_consensus_workflow(
        museq_vcf,
        samtools_vcf,
        rtg_vcf,
        freebayes_vcf,
        consensus_maf,
        chromosomes,
        sample_id,
        reference_vep
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='germline_consensus',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.germline_calling_consensus.consensus.main',
        args=(
            mgd.InputFile(museq_vcf),
            mgd.InputFile(freebayes_vcf),
            mgd.InputFile(rtg_vcf),
            mgd.InputFile(samtools_vcf),
            mgd.TempOutputFile('consensus.vcf'),
            mgd.TempOutputFile('counts.csv'),
            chromosomes
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
        kwargs={'normal_id': sample_id}
    )

    workflow.transform(
        name='maf_counts',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.germline_calling_consensus.tasks.update_maf_counts',
        args=(
            mgd.TempInputFile('consensus.maf'),
            mgd.TempInputFile('counts.csv'),
            mgd.OutputFile(consensus_maf),
        )
    )

    return workflow
