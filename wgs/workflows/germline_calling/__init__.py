import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_germline_calling_workflow(
        samples,
        normals,
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
        normal_ids,
        single_node=False
):
    museq_ss_vcf = dict([(sampid, museq_ss_vcf[sampid])
                         for sampid in samples])
    museq_ss_maf = dict([(sampid, museq_ss_maf[sampid])
                         for sampid in samples])
    museq_single_pdf = dict([(sampid, museq_single_pdf[sampid])
                             for sampid in samples])

    samtools_germline_vcf = dict([(sampid, samtools_germline_vcf[sampid])
                                  for sampid in samples])
    samtools_germline_maf = dict([(sampid, samtools_germline_maf[sampid])
                                  for sampid in samples])
    roh_calls = dict([(sampid, roh_calls[sampid])
                      for sampid in samples])

    freebayes_germline_vcf = dict([(sampid, freebayes_germline_vcf[sampid])
                                   for sampid in samples])
    freebayes_germline_maf = dict([(sampid, freebayes_germline_maf[sampid])
                                   for sampid in samples])

    rtg_germline_vcf = dict([(sampid, rtg_germline_vcf[sampid])
                             for sampid in samples])
    rtg_germline_maf = dict([(sampid, rtg_germline_maf[sampid])
                             for sampid in samples])

    consensus_germline_maf = dict([(sampid, consensus_germline_maf[sampid])
                                   for sampid in samples])

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

    workflow.subworkflow(
        name="mutationseq_single",
        func='wgs.workflows.mutationseq.create_museq_workflow',
        axes=('sample_id',),
        args=(
            mgd.OutputFile(
                'museq_germlines.vcf.gz', 'sample_id',
                extensions=['.csi', '.tbi'],
                fnames=museq_ss_vcf
            ),
            mgd.OutputFile(
                'museq_germlines.maf', 'sample_id',
                fnames=museq_ss_maf
            ),
            mgd.OutputFile('museq_single_pdf', 'sample_id', fnames=museq_single_pdf),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
        ),
        kwargs={
            'tumour_id': None,
            'normal_id': mgd.TempInputObj('normal_id', 'sample_id'),
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
            'single_node': single_node,
            'germline_refdata': paths_refdir['germline_portrait_ref'],
            'thousand_genomes': paths_refdir['thousand_genomes'],
            'dbsnp': paths_refdir['dbsnp'],
        }
    )

    workflow.subworkflow(
        name="samtools_germline",
        func='wgs.workflows.samtools_germline.create_samtools_germline_workflow',
        axes=('sample_id',),
        args=(
            mgd.OutputFile("samtools_germlines_anno.vcf.gz", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=samtools_germline_vcf),
            mgd.OutputFile("samtools_germlines_anno.maf", 'sample_id',
                           fnames=samtools_germline_maf),
            mgd.OutputFile("roh_calls.csv.gz", 'sample_id',
                           fnames=roh_calls, extensions=['.yaml']),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
            mgd.TempInputObj('normal_id', 'sample_id', fnames=normal_ids),
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="freebayes_germline",
        func='wgs.workflows.freebayes.create_freebayes_germline_workflow',
        axes=('sample_id',),
        args=(
            mgd.OutputFile("freebayes_germlines_anno.vcf.gz", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=freebayes_germline_vcf),
            mgd.OutputFile("freebayes_germlines_anno.maf", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=freebayes_germline_maf),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            paths_refdir['reference'],
            paths_refdir['reference_vep'],
            chromosomes,
            mgd.TempInputObj('normal_id', 'sample_id'),
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="rtg_germline",
        func='wgs.workflows.rtg_germline.create_rtg_germline_workflow',
        axes=('sample_id',),
        args=(
            mgd.OutputFile("rtg_germlines_anno.vcf.gz", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=rtg_germline_vcf),
            mgd.OutputFile("rtg_germlines_anno.maf", 'sample_id', extensions=['.csi', '.tbi'],
                           fnames=rtg_germline_maf),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            paths_refdir['reference'],
            paths_refdir['reference_sdf'],
            paths_refdir['reference_vep'],
            chromosomes,
            mgd.TempInputObj('normal_id', 'sample_id'),
        ),
        kwargs={
            'single_node': single_node,
        }
    )

    workflow.subworkflow(
        name="germline_consensus",
        func='wgs.workflows.germline_calling_consensus.create_germline_consensus_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('museq_germlines.vcf.gz', 'sample_id', fnames=museq_ss_vcf),
            mgd.InputFile("samtools_germlines_anno.vcf.gz", 'sample_id',
                          fnames=samtools_germline_vcf),
            mgd.InputFile("rtg_germlines_anno.vcf.gz", 'sample_id',
                          fnames=rtg_germline_vcf),
            mgd.InputFile("freebayes_germlines_anno.vcf.gz", 'sample_id',
                          fnames=freebayes_germline_vcf),
            mgd.OutputFile("germlines_consensus.maf", 'sample_id',
                           fnames=consensus_germline_maf),
            chromosomes,
            paths_refdir['reference_vep'],
            mgd.TempInputObj('normal_id', 'sample_id')
        ),
    )

    return workflow
