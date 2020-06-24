'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd

from wgs.config import config


def create_vcf2maf_workflow(
        museq_snv_vcf,
        strelka_snv_vcf,
        strelka_indel_vcf,
        maf_file,
        reference,
        tumour_id=None,
        normal_id=None
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    workflow.transform(
        name='overlap_museq_strelka_snv',
        func='wgs.workflows.vcf2maf.tasks.filter_vcfs',
        args=(
            mgd.InputFile(museq_snv_vcf),
            mgd.InputFile(strelka_snv_vcf),
            mgd.TempOutputFile('snvs.vcf'),
        )
    )

    workflow.transform(
        name='add_indels',
        func='wgs.utils.vcfutils.concatenate_vcf',
        args=(
            [mgd.TempInputFile('snvs.vcf'),
             mgd.InputFile(strelka_indel_vcf)],
            mgd.TempOutputFile('allcalls.vcf'),
        )
    )

    workflow.transform(
        name='vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        args=(
            mgd.TempInputFile('allcalls.vcf'),
            mgd.OutputFile(maf_file),
            mgd.TempSpace('vcf2maf_temp'),
            reference
        ),
        kwargs={
            'docker_image': config.containers('vcf2maf'),
            'tumour_id': tumour_id,
            'normal_id': normal_id
        }
    )

    return workflow
