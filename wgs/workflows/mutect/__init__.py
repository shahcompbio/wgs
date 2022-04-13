'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_mutect_workflow(
        normal_bam,
        tumour_bam,
        snv_vcf,
        snv_maf,
        reference,
        reference_vep,
        params_refdir,
        normal_id,
        tumour_id,
        single_node=None
):
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.mutect.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='1:00',
        ),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            params_refdir['chromosomes']
        ),
        kwargs={'size': params['split_size']}
    )

    if single_node:
        workflow.transform(
            name='mutect_one_node',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            func='wgs.workflows.mutect.tasks.run_mutect_one_job',
            args=(
                mgd.TempSpace("run_mutect_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(normal_bam),
                mgd.InputFile(tumour_bam)
            ),
        )
    else:
        workflow.transform(
            name='mutect_caller',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.mutect.tasks.run_mutect',
            args=(
                mgd.TempOutputFile('mutect.vcf', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(normal_bam),
                mgd.InputFile(tumour_bam),
                mgd.TempSpace('mutect_temp', 'interval')
            ),
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='8:00',
            ),
            func='wgs.workflows.mutect.tasks.merge_vcfs',
            args=(
                mgd.TempInputFile('mutect.vcf', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            ),
        )

    workflow.transform(
        name='bcftools_normalize',
        ctx=helpers.get_default_ctx(
            walltime='8:00',
        ),
        func='wgs.utils.vcfutils.bcftools_normalize',
        args=(
            mgd.TempInputFile('merged.vcf'),
            mgd.TempOutputFile('normalized.vcf'),
            reference,
        )
    )

    workflow.transform(
        name='finalise_snvs',
        ctx=helpers.get_default_ctx(
            walltime='8:00',
        ),
        func='wgs.utils.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('normalized.vcf'),
            mgd.OutputFile(snv_vcf, extensions=['.tbi', '.csi']),
        ),
    )

    workflow.subworkflow(
        name="strelka_indel_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(snv_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(snv_maf),
            reference_vep,
            params_refdir['vep_fasta_suffix'],
            params_refdir['ncbi_build'],
            params_refdir['cache_version']
        ),
        kwargs={'tumour_id': tumour_id, 'normal_id': normal_id}
    )

    return workflow
