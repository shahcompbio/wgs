'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_annotation_workflow(
        input_vcf,
        annotated_vcf,
        snpeff,
        mutationassessor,
        dbsnp,
        thousand_genomes,
        cosmic,
        mappability,
        chromosomes,
        vcftools_docker=None,
        snpeff_docker=None,
        input_type='snv',
):
    databases = {
        'snpeff_params': {'snpeff_config': snpeff, },
        'mutation_assessor_params': {'db': mutationassessor},
        'dbsnp_params': {'db': dbsnp},
        'thousandgen_params': {'db': thousand_genomes},
        'cosmic_params': {'db': cosmic},
        'mappability_ref': mappability
    }

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('chrom'),
        value=chromosomes,
    )

    workflow.transform(
        name='split_vcf_by_chrom',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.split_by_chrom',
        axes=('chrom',),
        args=(
            mgd.InputFile(input_vcf),
            mgd.TempOutputFile('split_chrom.vcf', 'chrom'),
            mgd.InputInstance('chrom')
        ),
    )

    workflow.transform(
        name='run_snpeff',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_snpeff',
        axes=('chrom',),
        args=(
            mgd.TempInputFile('split_chrom.vcf', 'chrom'),
            mgd.TempOutputFile('annotSnpEff.vcf', 'chrom'),
            databases,
        ),
        kwargs={'docker_image': snpeff_docker}
    )

    workflow.transform(
        name='run_mutation_assessor',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        axes=('chrom',),
        func='wgs.workflows.vcf_annotation.tasks.run_mutation_assessor',
        args=(
            mgd.TempInputFile('annotSnpEff.vcf', 'chrom'),
            mgd.TempOutputFile('annotMA.vcf', 'chrom'),
            databases,
            mgd.InputInstance('chrom')
        ),
    )

    workflow.transform(
        name='run_DBSNP',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        axes=('chrom',),
        func='wgs.workflows.vcf_annotation.tasks.run_DBSNP',
        args=(
            mgd.TempInputFile('annotMA.vcf', 'chrom'),
            mgd.TempOutputFile('flagDBsnp.vcf', 'chrom'),
            databases,
            mgd.InputInstance('chrom'),
        ),
        kwargs={'input_type': input_type},
    )

    workflow.transform(
        name='run_1000gen',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        axes=('chrom',),
        func='wgs.workflows.vcf_annotation.tasks.run_1000gen',
        args=(
            mgd.TempInputFile('flagDBsnp.vcf', 'chrom'),
            mgd.TempOutputFile('flag1000gen.vcf', 'chrom'),
            databases,
            mgd.InputInstance('chrom'),
        ),
        kwargs={'input_type': input_type},
    )

    workflow.transform(
        name='run_cosmic',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        axes=('chrom',),
        func='wgs.workflows.vcf_annotation.tasks.run_cosmic',
        args=(
            mgd.TempInputFile('flag1000gen.vcf', 'chrom'),
            mgd.TempOutputFile('cosmic.vcf', 'chrom'),
            databases,
            mgd.InputInstance('chrom'),
        ),
        kwargs={'input_type': input_type},
    )

    workflow.transform(
        name='low_mappability_flag',
        func='wgs.workflows.vcf_annotation.tasks.flag_low_mappability',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        axes=('chrom',),
        args=(
            mgd.TempInputFile('cosmic.vcf', 'chrom'),
            mgd.TempOutputFile('low_mapp.vcf', 'chrom'),
            databases['mappability_ref'],
            mgd.InputInstance('chrom')
        ),
    ),

    workflow.transform(
        name='merge_chroms',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='8:00', ),
        func='wgs.utils.vcfutils.concatenate_vcf',
        args=(
            mgd.TempInputFile('low_mapp.vcf', 'chrom'),
            mgd.TempOutputFile('low_mapp.vcf'),
        )
    )

    workflow.transform(
        name='finalize',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='8:00', ),
        func='wgs.utils.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('low_mapp.vcf'),
            mgd.OutputFile(annotated_vcf, extensions=['.csi', '.tbi']),
        ),
        kwargs={'docker_image': vcftools_docker}
    )

    return workflow
