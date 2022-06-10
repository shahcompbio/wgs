'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_freebayes_germline_workflow(
        germline_vcf,
        germline_maf,
        bam_file,
        reference,
        reference_vep,
        params_refdir,
        normal_id,
        single_node=None
):
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.freebayes.tasks.generate_intervals',
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
            name='freebayes_one_node',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            func='wgs.workflows.freebayes.tasks.run_freebayes_one_job',
            args=(
                mgd.TempSpace("run_freebayes_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(bam_file)
            )
        )
    else:
        workflow.transform(
            name='freebayes',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.freebayes.tasks.run_freebayes_germline',
            args=(
                mgd.TempOutputFile('freebayes_germline.vcf', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(bam_file),
                mgd.TempSpace('tempdir_freebayes', 'interval')
            ),
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='8:00',
            ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('freebayes_germline.vcf', 'interval'),
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
            mgd.OutputFile(germline_vcf, extensions=['.tbi', '.csi']),
        ),
    )

    workflow.subworkflow(
        name="freebayes_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_maf, extensions=['.tbi', '.csi']),
            reference_vep,
            params_refdir['vep_fasta_suffix'],
            params_refdir['ncbi_build'],
            params_refdir['cache_version'],
            params_refdir['species']
        ),
        kwargs={'normal_id': normal_id}
    )

    return workflow
