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
        params_refdir,
        reference_vep,
        normal_id
):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=params_refdir['chromosomes'],
    )

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
            params_refdir['chromosomes']
        ),
    )

    workflow.subworkflow(
        name="consensus_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.TempInputFile('consensus.vcf'),
            mgd.TempOutputFile('consensus.maf'),
            reference_vep,
            params_refdir['vep_fasta_suffix'],
            params_refdir['ncbi_build'],
            params_refdir['cache_version'],
            params_refdir['species'],
        ),
        kwargs={'normal_id': normal_id}
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
