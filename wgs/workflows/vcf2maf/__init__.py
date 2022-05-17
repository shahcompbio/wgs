'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


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
        name='split_vcf',
        func='wgs.workflows.vcf2maf.tasks.split_vcf',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        args=(
            mgd.InputFile(vcf_file),
            mgd.TempOutputFile('split.vcf', 'split')
        ),
        kwargs={'lines_per_file': 5000}
    )

    workflow.transform(
        name='vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        axes=('split',),
        args=(
            mgd.TempInputFile('split.vcf', 'split'),
            mgd.TempOutputFile('maf_file.maf', 'split'),
            mgd.TempSpace('vcf2maf_temp', 'split'),
            reference,
            vep_fasta_suffix,
            vep_ncbi_build,
            vep_cache_version,
        )
    )

    workflow.transform(
        name='merge_maf',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf2maf.tasks.merge_mafs',
        args=(
            mgd.TempInputFile('maf_file.maf', 'split'),
            mgd.TempOutputFile('maf_file_merged.maf')
        )
    )

    workflow.transform(
        name='update_ids',
        func='wgs.workflows.vcf2maf.tasks.update_ids',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        args=(
            mgd.TempInputFile('maf_file_merged.maf'),
            tumour_id,
            normal_id,
            mgd.OutputFile(maf_file),
        )
    )

    return workflow
