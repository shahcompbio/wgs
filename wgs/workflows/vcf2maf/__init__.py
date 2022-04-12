'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

from wgs.config import config


def create_vcf2maf_workflow(
        vcf_file,
        maf_file,
        reference,
        vep_fasta_suffix,
        vep_ncbi_build,
        vep_cache_version,
        tumour_id=None,
        normal_id=None
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        args=(
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('maf_file.maf'),
            mgd.TempSpace('vcf2maf_temp'),
            reference,
            vep_fasta_suffix,
            vep_ncbi_build,
            vep_cache_version,
        )
    )

    workflow.transform(
        name='update_ids',
        func='wgs.workflows.vcf2maf.tasks.update_ids',
        args=(
            mgd.TempInputFile('maf_file.maf'),
            tumour_id,
            normal_id,
            mgd.OutputFile(maf_file),
        )
    )

    return workflow
