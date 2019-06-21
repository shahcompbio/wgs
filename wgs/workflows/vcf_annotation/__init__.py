'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import vcf_tasks

import tasks
from wgs.utils import helpers

def create_annotation_workflow(
        input_vcf,
        annotated_vcf,
        global_config,
        varcall_config,
        vcftools_docker=None,
        snpeff_docker=None
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='run_snpeff',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.run_snpeff,
        args=(
            mgd.InputFile(input_vcf),
            mgd.TempOutputFile('museq_annotSnpEff.vcf'),
            varcall_config,
        ),
        kwargs={'docker_image': snpeff_docker}
    )

    workflow.transform(
        name='run_mutation_assessor',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='8:00', ),
        func=tasks.run_mutation_assessor,
        args=(
            mgd.TempInputFile('museq_annotSnpEff.vcf'),
            mgd.TempOutputFile('museq_annotMA.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_DBSNP',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.run_DBSNP,
        args=(
            mgd.TempInputFile('museq_annotMA.vcf'),
            mgd.TempOutputFile('museq_flagDBsnp.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_1000gen',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.run_1000gen,
        args=(
            mgd.TempInputFile('museq_flagDBsnp.vcf'),
            mgd.TempOutputFile('museq_flag1000gen.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_cosmic',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=tasks.run_cosmic,
        args=(
            mgd.TempInputFile('museq_flag1000gen.vcf'),
            mgd.TempOutputFile('museq_cosmic.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='finalize',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func=vcf_tasks.finalise_vcf,
        args=(
            mgd.TempInputFile('museq_cosmic.vcf'),
            mgd.OutputFile(annotated_vcf),
        ),
        kwargs={'docker_image': vcftools_docker}
    )

    return workflow
