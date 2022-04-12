'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_samtools_germline_workflow(
        germline_vcf,
        germline_maf,
        germline_roh,
        bam_file,
        reference,
        reference_vep,
        chromosomes,
        normal_id,
        single_node=None
):
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.samtools_germline.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='1:00',
        ),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            chromosomes
        ),
        kwargs={'size': params['split_size']}
    )

    if single_node:
        workflow.transform(
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline_one_job',
            args=(
                mgd.TempSpace("run_samtools_temp"),
                mgd.TempOutputFile('merged.vcf'),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(bam_file)
            ),
        )
    else:
        workflow.transform(
            name='samtools_germline',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.samtools_germline.tasks.run_samtools_germline',
            args=(
                mgd.TempOutputFile('germline.vcf.gz', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(bam_file),
                mgd.TempSpace('tempdir_samtools', 'interval')
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
                mgd.TempInputFile('germline.vcf.gz', 'interval'),
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

    workflow.transform(
        name='roh_calling',
        ctx=helpers.get_default_ctx(
            walltime='8:00',
        ),
        func='wgs.workflows.samtools_germline.tasks.roh_calling',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_roh, extensions=['.yaml']),
            mgd.TempSpace('roh_calling_temp')
        ),
    )

    workflow.subworkflow(
        name="samtools_germline_maf",
        func='wgs.workflows.vcf2maf.create_vcf2maf_workflow',
        args=(
            mgd.InputFile(germline_vcf, extensions=['.tbi', '.csi']),
            mgd.OutputFile(germline_maf, extensions=['.tbi', '.csi']),
            reference_vep,
            params['vep_fasta_suffix'],
            params['ncbi_build'],
            params['cache_version']
        ),
        kwargs={'normal_id': normal_id}
    )

    return workflow
