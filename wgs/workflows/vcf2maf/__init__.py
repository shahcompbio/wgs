'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

from wgs.config import config
from wgs.utils import helpers

def create_vcf2maf_workflow(
        museq_snv_vcf,
        strelka_snv_vcf,
        strelka_indel_vcf,
        museq_germlines,
        samtools_germlines,
        maf_file,
        germline_maf_file,
        reference,
        chromosomes,
        sample_id,
        tumour_id=None,
        normal_id=None
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    workflow.transform(
        name='snv_consensus',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf2maf.consensus.main',
        args=(
            mgd.InputFile(museq_snv_vcf),
            mgd.InputFile(strelka_snv_vcf),
            mgd.InputFile(strelka_indel_vcf),
            mgd.TempOutputFile('consensus.vcf'),
            mgd.TempOutputFile('counts.csv'),
            chromosomes,
        ),
    )

    workflow.transform(
        name='vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        args=(
            mgd.TempInputFile('consensus.vcf'),
            mgd.TempOutputFile('consensus.maf'),
            mgd.TempSpace('vcf2maf_temp'),
            reference
        ),
        kwargs={
            'vcf2maf_docker_image': config.containers('vcf2maf'),
            'vcftools_docker_image': config.containers('vcftools'),
        }
    )

    workflow.transform(
        name='maf_counts',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf2maf.tasks.update_maf_counts',
        args=(
            mgd.TempInputFile('consensus.maf'),
            mgd.TempInputFile('counts.csv'),
            mgd.OutputFile(maf_file),
            sample_id,
        )
    )


    workflow.transform(
        name='germline_consensus',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf2maf.germline_consensus.main',
        args=(
            mgd.InputFile(museq_germlines),
            mgd.InputFile(samtools_germlines),
            mgd.TempOutputFile('germline_consensus.vcf'),
            mgd.TempOutputFile('germline_counts.csv'),
            chromosomes,
        ),
    )


    workflow.transform(
        name='germline_vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='48:00', ),
        args=(
            mgd.TempInputFile('germline_consensus.vcf'),
            mgd.TempOutputFile('germline_consensus.maf'),
            mgd.TempSpace('germline_vcf2maf_temp'),
            reference
        ),
        kwargs={
            'docker_image': config.containers('vcf2maf'),
            'tumour_id': tumour_id,
            'normal_id': normal_id
        }
    )

    workflow.transform(
        name='germline_maf_counts',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf2maf.tasks.update_maf_counts',
        args=(
            mgd.TempInputFile('germline_consensus.maf'),
            mgd.TempInputFile('germline_counts.csv'),
            mgd.OutputFile(germline_maf_file),
            sample_id,
        )
    )

    return workflow
