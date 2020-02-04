import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_titan_workflow(
        tumour_bam, normal_bam, targets, titan_raw_dir,
        segments, params, markers,
        global_config, config,
        intervals, sample_id,
        single_node=None
):
    museq_vcf = os.path.join(titan_raw_dir, '{}_museq.vcf'.format(sample_id))

    # optimal
    optimal_outfile = os.path.join(titan_raw_dir, '{}_titan_markers.csv.gz'.format(sample_id))
    optimal_params = os.path.join(titan_raw_dir, '{}_titan_params.csv.gz'.format(sample_id))
    optimal_segs = os.path.join(titan_raw_dir, '{}_titan_segs.csv.gz'.format(sample_id))
    optimal_igv_segs = os.path.join(titan_raw_dir, '{}_titan_igv_segs.csv.gz'.format(sample_id))
    optimal_parsed = os.path.join(titan_raw_dir, '{}_titan_parsed.csv.gz'.format(sample_id))
    optimal_plots = os.path.join(titan_raw_dir, '{}_titan_plots.pdf'.format(sample_id))

    tar_outputs = os.path.join(titan_raw_dir, '{}_data_all_parameters.tar.gz'.format(sample_id))

    chunks = [(v['num_clusters'], v['ploidy']) for v in intervals]

    targets = mgd.InputFile(targets) if targets else None

    ctx = {'docker_image': config['docker']['wgs']}

    workflow = pypeliner.workflow.Workflow(ctx=ctx)

    workflow.setobj(
        obj=mgd.OutputChunks('numclusters', 'ploidy'),
        value=chunks,
    )

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.titan.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='2:00', ),
        ret=mgd.OutputChunks('interval'),
        args=(
            config['reference_genome'],
            config['chromosomes']
        ),
        kwargs={'size': config['split_size']}
    )

    if single_node:
        workflow.transform(
            name='run_museq',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='96:00',
                ncpus=global_config['threads']),
            func='wgs.utils.museq_utils.run_museq_one_job',
            args=(
                mgd.TempSpace("run_museq_temp"),
                mgd.OutputFile(museq_vcf),
                config['reference_genome'],
                mgd.InputChunks('interval'),
                config['museq_params'],
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam, extensions=['.bai']),
                'normal_bam': mgd.InputFile(normal_bam, extensions=['.bai']),
                'titan_mode': True,
                'museq_docker_image': config['docker']['mutationseq'],
                'vcftools_docker_image': config['docker']['vcftools']
            }
        )
    else:
        workflow.transform(
            name='run_museq',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='24:00'),
            axes=('interval',),
            func='wgs.utils.museq_utils.run_museq',
            args=(
                mgd.TempOutputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('museq.log', 'interval'),
                config['reference_genome'],
                mgd.InputInstance('interval'),
                config['museq_params']
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam, extensions=['.bai']),
                'normal_bam': mgd.InputFile(normal_bam, extensions=['.bai']),
                'titan_mode': True,
                'docker_image': config['docker']['mutationseq']
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=global_config['memory']['high'],
                walltime='4:00', ),
            func='wgs.utils.museq_utils.merge_vcfs',
            args=(
                mgd.TempInputFile('museq.vcf', 'interval'),
                mgd.OutputFile(museq_vcf),
                mgd.TempSpace('merge_vcf'),
            ),
            kwargs={'docker_image': config['docker']['vcftools']}
        )

    workflow.transform(
        name='convert_museq_vcf2counts',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.convert_museq_vcf2counts',
        args=(
            mgd.InputFile(museq_vcf),
            mgd.TempOutputFile('museq_postprocess.txt'),
            config
        ),
    )

    workflow.transform(
        name='run_readcounter_tumour',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='16:00',
            disk=200
        ),
        func='wgs.workflows.titan.tasks.run_readcounter',
        args=(
            mgd.InputFile(tumour_bam, extensions=['.bai']),
            mgd.TempOutputFile('tumour.wig'),
            config,
        ),
    )

    workflow.transform(
        name='run_readcounter_normal',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='16:00',
            disk=200
        ),
        func='wgs.workflows.titan.tasks.run_readcounter',
        args=(
            mgd.InputFile(normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('normal.wig'),
            config,
        ),
    )

    workflow.transform(
        name='calc_correctreads_wig',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.calc_correctreads_wig',
        args=(
            mgd.TempInputFile('tumour.wig'),
            mgd.TempInputFile('normal.wig'),
            targets,
            mgd.TempOutputFile('correct_reads.txt'),
            config,
        ),
        kwargs={'docker_image': config['docker']['titan']}
    )

    workflow.transform(
        name='run_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='72:00',
            ncpus=config['ncpus']),
        func='wgs.workflows.titan.tasks.run_titan',
        args=(
            mgd.TempInputFile('museq_postprocess.txt'),
            mgd.TempInputFile('correct_reads.txt'),
            mgd.TempOutputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan.Rdata', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_params', 'numclusters', 'ploidy'),
            config['titan_params'],
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy'),
            sample_id
        ),
        kwargs={'docker_image': config['docker']['titan']}
    )

    workflow.transform(
        name='plot_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
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
            'chromosomes': config['chromosomes'],
            'docker_image': config['docker']['titan'],
        },
    )

    workflow.transform(
        name='calc_cnsegments_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.calc_cnsegments_titan',
        args=(
            mgd.TempInputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_igv', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('segs.csv', 'numclusters', 'ploidy'),
        ),
        kwargs={'docker_image': config['docker']['titan']}
    )

    workflow.transform(
        name='annot_pygenes',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.annot_pygenes',
        args=(
            mgd.TempInputFile('segs.csv', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_segs.csv', 'numclusters', 'ploidy'),
            config,
        ),
    )

    workflow.transform(
        name='parse_titan',
        axes=('numclusters', 'ploidy'),
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='4:00', ),
        func='wgs.workflows.titan.tasks.parse_titan_data',
        args=(
            mgd.TempInputFile('titan_segs.csv', 'numclusters', 'ploidy'),
            mgd.TempInputFile('titan_outfile', 'numclusters', 'ploidy'),
            mgd.TempOutputFile('titan_parsed.csv', 'numclusters', 'ploidy'),
            config['parse_titan'],
        ),
    )

    # select optimal solution
    workflow.transform(
        name="select_optimal_solution",
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
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
            mgd.OutputFile(optimal_segs, extensions=['.yaml']),
            mgd.OutputFile(optimal_igv_segs, extensions=['.yaml']),
            mgd.OutputFile(optimal_params, extensions=['.yaml']),
            mgd.OutputFile(optimal_outfile, extensions=['.yaml']),
            mgd.OutputFile(optimal_parsed, extensions=['.yaml']),
            mgd.OutputFile(optimal_plots),
        )
    )

    workflow.transform(
        name='tar_all_data',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
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
