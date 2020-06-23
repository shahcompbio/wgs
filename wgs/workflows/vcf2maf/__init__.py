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
        tumour_id=None,
        normal_id=None
):
    workflow = pypeliner.workflow.Workflow(
        ctx={'docker_image': config.containers('wgs')}
    )

    workflow.transform(
        name='vcf2maf',
        func='wgs.workflows.vcf2maf.tasks.run_vcf2maf',
        args=(
            mgd.InputFile(vcf_file),
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
