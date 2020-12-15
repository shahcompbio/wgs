import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_somatic_calling_workflow(
        samples,
        tumours,
        normals,
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
        normal_ids,
        tumour_ids,
        single_node=False,
        is_exome=False
):
    strelka_snv_vcf = dict([(sampid, strelka_snv_vcf[sampid])
                            for sampid in samples])
    strelka_indel_vcf = dict([(sampid, strelka_indel_vcf[sampid])
                              for sampid in samples])
    strelka_snv_maf = dict([(sampid, strelka_snv_maf[sampid])
                            for sampid in samples])
    strelka_indel_maf = dict([(sampid, strelka_indel_maf[sampid])
                              for sampid in samples])

    museq_vcf = dict([(sampid, museq_vcf[sampid])
                      for sampid in samples])
    museq_maf = dict([(sampid, museq_maf[sampid])
                      for sampid in samples])
    museq_paired_pdf = dict([(sampid, museq_paired_pdf[sampid])
                             for sampid in samples])

    mutect_vcf = dict([(sampid, mutect_vcf[sampid])
                       for sampid in samples])
    mutect_maf = dict([(sampid, mutect_maf[sampid])
                       for sampid in samples])

    somatic_consensus_maf = dict([(sampid, somatic_consensus_maf[sampid]) for sampid in samples])

    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    paths_refdir = config.refdir_data(refdir)['paths']

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config.containers('wgs'))
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.setobj(
        obj=mgd.TempOutputObj('normal_id', 'sample_id', axes_origin=[]),
        value={v: normal_ids[v] for v in samples})

    workflow.setobj(
        obj=mgd.TempOutputObj('tumour_id', 'sample_id', axes_origin=[]),
        value={v: tumour_ids[v] for v in samples})

    workflow.subworkflow(
        name="mutationseq_paired",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        axes=('sample_id',),
        args=(
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=museq_vcf),
            mgd.OutputFile('museq_snv_ann.maf', 'sample_id',
                           fnames=museq_maf),
            mgd.OutputFile('museq_paired_pdf', 'sample_id', fnames=museq_paired_pdf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
        ),
        kwargs={
            'normal_id': mgd.TempInputObj('normal_id', 'sample_id'),
            'tumour_id': mgd.TempInputObj('tumour_id', 'sample_id'),
            'tumour_bam': mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                                        extensions=['.bai'], axes_origin=[]),
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="strelka",
        func='wgs.workflows.strelka.create_strelka_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.OutputFile('strelka_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_snv_vcf),
            mgd.OutputFile('strelka_snv_ann.maf', 'sample_id',
                           fnames=strelka_snv_maf),
            mgd.OutputFile('strelka_indel_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=strelka_indel_vcf),
            mgd.OutputFile('strelka_indel_ann.maf', 'sample_id',
                           fnames=strelka_indel_maf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
            mgd.TempInputObj('normal_id', 'sample_id'),
            mgd.TempInputObj('tumour_id', 'sample_id'),
        ),
        kwargs={
            'single_node': single_node,
            'is_exome': is_exome
        },

    )

    workflow.subworkflow(
        name="mutect",
        func='wgs.workflows.mutect.create_mutect_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.OutputFile('mutect_snv_ann.vcf.gz', 'sample_id',
                           extensions=['.csi', '.tbi'], fnames=mutect_vcf),
            mgd.OutputFile('mutect_snv_ann.maf', 'sample_id',
                           fnames=mutect_maf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
            mgd.TempInputObj('normal_id', 'sample_id'),
            mgd.TempInputObj('tumour_id', 'sample_id'),
        ),
        kwargs={
            'single_node': single_node,
        },

    )

    workflow.subworkflow(
        name="somatic_consensus",
        func='wgs.workflows.somatic_calling_consensus.create_somatic_consensus_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('mutect_snv_ann.vcf.gz', 'sample_id',
                          extensions=['.csi', '.tbi'], fnames=mutect_vcf),
            mgd.InputFile('strelka_snv_ann.vcf.gz', 'sample_id',
                          extensions=['.csi', '.tbi'], fnames=strelka_snv_vcf),
            mgd.InputFile('strelka_indel_ann.vcf.gz', 'sample_id',
                          extensions=['.csi', '.tbi'], fnames=strelka_indel_vcf),
            mgd.InputFile('museq_snv_ann.vcf.gz', 'sample_id',
                          extensions=['.csi', '.tbi'], fnames=museq_vcf),
            mgd.OutputFile("somatic_consensus.maf", 'sample_id',
                           fnames=somatic_consensus_maf),
            chromosomes,
            paths_refdir['reference_vep'],
            mgd.TempInputObj('normal_id', 'sample_id'),
            mgd.TempInputObj('tumour_id', 'sample_id'),
        ),
    )

    return workflow
