import pypeliner
import pypeliner.managed as mgd
from wgs.config import config


def create_somatic_calling_workflow(
        tumour,
        normal,
        museq_vcf,
        museq_maf,
        museq_paired_pdf,
        strelka_snv_vcf,
        strelka_snv_maf,
        strelka_indel_vcf,
        strelka_indel_maf,
        mutect_vcf,
        mutect_maf,
        somatic_consensus_maf,
        refdir,
        normal_id,
        tumour_id,
        single_node=False,
        is_exome=False
):
    paths_refdir = config.refdir_data(refdir)['paths']
    params_refdir = config.refdir_data(refdir)['params']

    workflow = pypeliner.workflow.Workflow()

    workflow.subworkflow(
        name="mutationseq_paired",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        args=(
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_vcf),
            mgd.OutputFile('museq_snv_ann.maf', 'sample_id',
                           fnames=museq_maf),
            mgd.OutputFile('museq_paired_pdf', 'sample_id', fnames=museq_paired_pdf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
        ),
        kwargs={
            'normal_id': normal_id,
            'tumour_id': tumour_id,
            'tumour_bam': mgd.InputFile(tumour, extensions=['.bai']),
            'normal_bam': mgd.InputFile(normal, extensions=['.bai']),
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="strelka",
        func='wgs.workflows.strelka.create_strelka_workflow',
        args=(
            mgd.InputFile(normal, extensions=['.bai']),
            mgd.InputFile(tumour, extensions=['.bai']),
            mgd.OutputFile(strelka_snv_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(strelka_snv_maf),
            mgd.OutputFile(strelka_indel_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(strelka_indel_maf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
            normal_id,
            tumour_id,
        ),
        kwargs={
            'single_node': single_node,
            'is_exome': is_exome
        },
    )

    workflow.subworkflow(
        name="mutect",
        func='wgs.workflows.mutect.create_mutect_workflow',
        args=(
            mgd.InputFile(normal, extensions=['.bai']),
            mgd.InputFile(tumour, extensions=['.bai']),
            mgd.OutputFile(mutect_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(mutect_maf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            params_refdir,
            normal_id,
            tumour_id,
        ),
        kwargs={
            'single_node': single_node,
        },

    )

    workflow.subworkflow(
        name="somatic_consensus",
        func='wgs.workflows.somatic_calling_consensus.create_somatic_consensus_workflow',
        args=(
            mgd.InputFile(mutect_vcf, extensions=['.csi', '.tbi']),
            mgd.InputFile(strelka_snv_vcf, extensions=['.csi', '.tbi']),
            mgd.InputFile(strelka_indel_vcf, extensions=['.csi', '.tbi']),
            mgd.InputFile(museq_vcf, extensions=['.csi', '.tbi']),
            mgd.OutputFile(somatic_consensus_maf),
            params_refdir,
            paths_refdir['reference_vep'],
            normal_id,
            tumour_id,
        ),
    )

    return workflow
