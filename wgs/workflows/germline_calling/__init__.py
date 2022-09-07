import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


def create_germline_calling_workflow(
        normal,
        museq_ss_vcf,
        museq_ss_maf,
        museq_single_pdf,
        samtools_germline_vcf,
        samtools_germline_maf,
        roh_calls,
        freebayes_germline_vcf,
        freebayes_germline_maf,
        rtg_germline_vcf,
        rtg_germline_maf,
        consensus_germline_maf,
        refdir,
        normal_id,
        single_node=False
):
    params_refdir = config.refdir_data(refdir)['params']
    paths_refdir = config.refdir_data(refdir)['paths']

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name="mutationseq_single",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        args=(
            mgd.OutputFile(museq_ss_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(museq_ss_maf),
            mgd.OutputFile(museq_single_pdf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
        ),
        kwargs={
            'tumour_id': None,
            'normal_id': normal_id,
            'tumour_bam': None,
            'normal_bam': mgd.InputFile(normal, extensions=['.bai']),
            'single_node': single_node,
            'germline_refdata': paths_refdir['germline_portrait_ref'],
            'thousand_genomes': paths_refdir['thousand_genomes'],
            'dbsnp': paths_refdir['dbsnp'],
        }
    )

    workflow.subworkflow(
        name="samtools_germline",
        func='wgs.workflows.samtools_germline.create_samtools_germline_workflow',
        args=(
            mgd.OutputFile(samtools_germline_vcf),
            mgd.OutputFile(samtools_germline_maf),
            mgd.OutputFile(roh_calls, extensions=['.yaml']),
            mgd.InputFile(normal, extensions=['.bai']),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
            normal_id,
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="freebayes_germline",
        func='wgs.workflows.freebayes.create_freebayes_germline_workflow',
        args=(
            mgd.OutputFile(freebayes_germline_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(freebayes_germline_maf, extensions=['.csi', '.tbi']),
            mgd.InputFile(normal, extensions=['.bai']),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
            normal_id,
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="rtg_germline",
        func='wgs.workflows.rtg_germline.create_rtg_germline_workflow',
        args=(
            mgd.OutputFile(rtg_germline_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(rtg_germline_maf, extensions=['.csi', '.tbi'], ),
            mgd.InputFile(normal, extensions=['.bai']),
            paths_refdir['reference'],
            paths_refdir['reference_sdf'],
            paths_refdir['reference_vep'],
            params_refdir,
            normal_id,
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="germline_consensus",
        func='wgs.workflows.germline_calling_consensus.create_germline_consensus_workflow',
        args=(
            mgd.InputFile(museq_ss_vcf),
            mgd.InputFile(samtools_germline_vcf),
            mgd.InputFile(rtg_germline_vcf),
            mgd.InputFile(freebayes_germline_vcf),
            mgd.OutputFile(consensus_germline_maf),
            params_refdir,
            paths_refdir['reference_vep'],
            normal_id
        ),
    )

    return workflow
