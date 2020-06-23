import os

import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def create_titan_workflow(
        tumour_bam, normal_bam, targets, outfile, params, segs, igv_segs,
        parsed, plots, tar_outputs, museq_vcf,
        sample_id, reference, chromosomes, het_positions, map_wig, gc_wig, pygenes_gtf,
        single_node=None
):
    cn_params = config.default_params('copynumber_calling')

    chunks = [(v['num_clusters'], v['ploidy']) for v in cn_params['titan_intervals']]

    targets = mgd.InputFile(targets) if targets else None

    ctx = {'docker_image': config.containers('wgs')}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('numclusters', 'ploidy'),
        value=chunks,
    )

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.titan.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='2:00', ),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            chromosomes,
        ),
        kwargs={'size': cn_params['split_size']}
    )

    if single_node:
        workflow.transform(
            name='run_museq',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='96:00',
                ncpus=8),
            func='wgs.utils.museq_utils.run_museq_one_job',
            args=(
                mgd.TempSpace("run_museq_temp"),
                mgd.OutputFile(museq_vcf),
                reference,
                mgd.InputChunks('interval'),
                cn_params['museq_params'],
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam, extensions=['.bai']),
                'normal_bam': mgd.InputFile(normal_bam, extensions=['.bai']),
                'titan_mode': True,
                'museq_docker_image': config.containers('mutationseq'),
                'vcftools_docker_image': config.containers('vcftools')
            }
        )
    else:
        workflow.transform(
            name='run_museq',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00'),
            axes=('interval',),
            func='wgs.utils.museq_utils.run_museq',
            args=(
                mgd.TempOutputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('museq.log', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                cn_params['museq_params'],
                mgd.TempSpace('museq_temp_titan', 'interval')
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam, extensions=['.bai']),
                'normal_bam': mgd.InputFile(normal_bam, extensions=['.bai']),
                'titan_mode': True,
                'docker_image': config.containers('mutationseq')
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='4:00', ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('museq.vcf', 'interval'),
                mgd.OutputFile(museq_vcf),
                mgd.TempSpace('merge_vcf'),
            ),
            kwargs={'docker_image': config.containers('vcftools')}
        )

    workflow.transform(
        name='convert_museq_vcf2counts',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.convert_museq_vcf2counts',
        args=(
            mgd.InputFile(museq_vcf),
            mgd.TempOutputFile('museq_postprocess.txt'),
            het_positions,
        ),
    )

    workflow.transform(
        name='run_readcounter_tumour',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='16:00',
            disk=200
        ),
        func='wgs.workflows.titan.tasks.run_readcounter',
        args=(
            mgd.InputFile(tumour_bam, extensions=['.bai']),
            mgd.TempOutputFile('tumour.wig'),
            chromosomes,
            cn_params['readcounter']
        ),
    )

    workflow.transform(
        name='run_readcounter_normal',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='16:00',
            disk=200
        ),
        func='wgs.workflows.titan.tasks.run_readcounter',
        args=(
            mgd.InputFile(normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('normal.wig'),
            chromosomes,
            cn_params['readcounter']
        ),
    )

    workflow.transform(
        name='calc_correctreads_wig',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.calc_correctreads_wig',
        args=(
            mgd.TempInputFile('tumour.wig'),
            mgd.TempInputFile('normal.wig'),
            targets,
            mgd.TempOutputFile('correct_reads.txt'),
            gc_wig,
            map_wig,
            cn_params['genome_type']
        ),
        kwargs={'docker_image': config.containers('titan')}
    )

    workflow.transform(
        name='run_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='72:00',
            ncpus='8'),
        func='wgs.workflows.titan.tasks.run_titan',
        args=(
            mgd.TempInputFile('museq_postprocess.txt'),
            mgd.TempInputFile('correct_reads.txt'),
            mgd.TempOutputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan.Rdata', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_params', 'numclusters', 'ploidy'),
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy'),
            sample_id,
            map_wig,
            cn_params['titan_params'],
            cn_params['genome_type']
        ),
        kwargs={'docker_image': config.containers('titan'), 'threads': '8'}
    )

    workflow.transform(
        name='plot_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='16:00', ),
        func='wgs.workflows.titan.tasks.plot_titan',
        args=(
            mgd.TempInputFile('titan.Rdata', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_plots', 'numclusters', 'ploidy'),
            mgd.TempSpace("titan_plots_tempdir", 'numclusters', 'ploidy'),
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy')
        ),
        kwargs={
            'chromosomes': chromosomes,
            'docker_image': config.containers('titan'),
        },
    )

    workflow.transform(
        name='calc_cnsegments_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.calc_cnsegments_titan',
        args=(
            mgd.TempInputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_igv', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('segs.csv', 'numclusters', 'ploidy'),
            sample_id,
        ),
        kwargs={'docker_image': config.containers('titan')}
    )

    workflow.transform(
        name='annot_pygenes',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.annot_pygenes',
        args=(
            mgd.TempInputFile('segs.csv', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_segs.csv', 'numclusters', 'ploidy'),
            pygenes_gtf,
        ),
    )

    workflow.transform(
        name='parse_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.parse_titan_data',
        args=(
            mgd.TempInputFile('titan_segs.csv', 'numclusters', 'ploidy'),
            mgd.TempInputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_parsed.csv', 'numclusters', 'ploidy'),
        ),
    )

    # select optimal solution
    workflow.transform(
        name="select_optimal_solution",
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00', ),
        func="wgs.workflows.titan.tasks.select_optimal_solution",
        args=(
            chunks,
            mgd.TempInputFile('titan_params', 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_segs.csv", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile('titan_igv', 'numclusters', 'ploidy'),
            mgd.TempInputFile("titan_outfile", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_parsed.csv", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_plots", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.OutputFile(segs, extensions=['.yaml']),
            mgd.OutputFile(igv_segs, extensions=['.yaml']),
            mgd.OutputFile(params, extensions=['.yaml']),
            mgd.OutputFile(outfile, extensions=['.yaml']),
            mgd.OutputFile(parsed, extensions=['.yaml']),
            mgd.OutputFile(plots),
        )
    )

    workflow.transform(
        name='tar_all_data',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='4:00', ),
        func="wgs.workflows.titan.tasks.tar_all_data",
        args=(
            mgd.TempInputFile('titan_params', 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_segs.csv", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile('titan_igv', 'numclusters', 'ploidy'),
            mgd.TempInputFile("titan_outfile", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_parsed.csv", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.TempInputFile("titan_plots", 'numclusters', 'ploidy', axes_origin=[]),
            mgd.OutputFile(tar_outputs),
            mgd.TempSpace("titan_all_parameters_data"),
            chunks
        )
    )

    return workflow
