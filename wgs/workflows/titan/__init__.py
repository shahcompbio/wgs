import os

import pypeliner
import pypeliner.managed as mgd

import tasks
import wgs.utils.museq_utils

def create_titan_workflow(
        tumour_bam, normal_bam, targets, titan_raw_dir,
        segments, params, markers,
        global_config, config,
        intervals, sample_id,
        single_node=None
):
    titan_outdir = os.path.join(titan_raw_dir, 'clusters_{numclusters}', 'ploidy_{ploidy}')
    igv_template = os.path.join(titan_outdir, 'igv_segs.txt')
    outfile_template = os.path.join(titan_outdir, 'titan_markers.txt')
    params_template = os.path.join(titan_outdir, 'titan_params.txt')
    segs_template = os.path.join(titan_outdir, 'titan_segs.txt')
    plots_template = os.path.join(titan_outdir, 'titan_plots.tar.gz')
    parsed_template = os.path.join(titan_outdir, 'titan_parsed.csv')
    museq_vcf = os.path.join(titan_raw_dir, 'museq.vcf')

    chunks = [(v['num_clusters'], v['ploidy']) for v in intervals]

    targets = mgd.InputFile(targets) if targets else None

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('numclusters', 'ploidy'),
        value=chunks,
    )

    workflow.transform(
        name='generate_intervals',
        func=tasks.generate_intervals,
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '02:00'},
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
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '24:00'},
            func=wgs.utils.museq_utils.run_museq_one_job,
            args=(
                mgd.TempSpace("run_museq_temp"),
                mgd.TempOutputFile('merged.vcf'),
                config['reference_genome'],
                mgd.InputChunks('interval'),
                config['museq_params'],
            ),
            kwargs={
                'tumour_bam': mgd.InputFile(tumour_bam, extensions=['.bai']),
                'normal_bam': mgd.InputFile(normal_bam, extensions=['.bai']),
                'titan_mode': True
            }
        )
    else:
        workflow.transform(
            name='run_museq',
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '24:00'},
            axes=('interval',),
            func=wgs.utils.museq_utils.run_museq,
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
                'titan_mode': True
            }
        )

        workflow.transform(
            name='merge_vcfs',
            ctx={'num_retry': 3, 'mem_retry_increment': 2,
                 'mem': global_config['memory']['high'],
                 'ncpus': global_config['threads'],
                 'walltime': '08:00'},
            func=tasks.merge_vcfs,
            args=(
                mgd.TempInputFile('museq.vcf', 'interval'),
                mgd.TempOutputFile('merged.vcf'),
                mgd.TempSpace('merge_vcf'),
            )
        )

    workflow.transform(
        name='convert_museq_vcf2counts',
        ctx={'mem': global_config['memory']['high'],
             'ncpus': 1, 'walltime': '02:00'},
        func=tasks.convert_museq_vcf2counts,
        args=(
            mgd.InputFile(museq_vcf),
            mgd.TempOutputFile('museq_postprocess.txt'),
            config,
        ),
    )

    workflow.transform(
        name='run_readcounter_tumour',
        ctx={'mem': global_config['memory']['high'],
             'ncpus': 1, 'walltime': '02:00'},
        func=tasks.run_readcounter,
        args=(
            mgd.InputFile(tumour_bam, extensions=['.bai']),
            mgd.TempOutputFile('tumour.wig'),
            config,
        ),
    )

    workflow.transform(
        name='run_readcounter_normal',
        ctx={'mem': global_config['memory']['high'],
             'ncpus': 1, 'walltime': '02:00'},
        func=tasks.run_readcounter,
        args=(
            mgd.InputFile(normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('normal.wig'),
            config,
        ),
    )

    workflow.transform(
        name='calc_correctreads_wig',
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '02:00'},
        func=tasks.calc_correctreads_wig,
        args=(
            mgd.TempInputFile('tumour.wig'),
            mgd.TempInputFile('normal.wig'),
            targets,
            mgd.TempOutputFile('correct_reads.txt'),
            config,
        ),
    )

    workflow.transform(
        name='run_titan',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['high'],
             'ncpus': 1, 'walltime': '06:00'},
        func=tasks.run_titan,
        args=(
            mgd.TempInputFile('museq_postprocess.txt'),
            mgd.TempInputFile('correct_reads.txt'),
            mgd.OutputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
            mgd.TempOutputFile('titan.Rdata', 'numclusters', 'ploidy'),
            mgd.OutputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
            config['titan_params'],
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy')
        )
    )

    workflow.transform(
        name='plot_titan',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['low'],
             'ncpus': 1, 'walltime': '02:00'},
        func=tasks.plot_titan,
        args=(
            mgd.TempInputFile('titan.Rdata', 'numclusters', 'ploidy'),
            mgd.InputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
            mgd.OutputFile('titan_plots', 'numclusters', 'ploidy', template=plots_template),
            mgd.TempSpace("titan_plots_tempdir", 'numclusters', 'ploidy'),
            config,
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy')
        ),
    )
    #
    # workflow.transform(
    #     name='calc_cnsegments_titan',
    #     axes=('numclusters', 'ploidy'),
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.calc_cnsegments_titan,
    #     args=(
    #         mgd.InputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
    #         mgd.OutputFile('titan_igv', 'numclusters', 'ploidy', template=igv_template),
    #         mgd.TempOutputFile('segs.csv', 'numclusters', 'ploidy'),
    #     ),
    # )
    #
    # workflow.transform(
    #     name='annot_pygenes',
    #     axes=('numclusters', 'ploidy'),
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.annot_pygenes,
    #     args=(
    #         mgd.TempInputFile('segs.csv', 'numclusters', 'ploidy'),
    #         mgd.OutputFile('titan_segs.csv', 'numclusters', 'ploidy', template=segs_template),
    #         config,
    #     ),
    # )
    #
    # workflow.transform(
    #     name='parse_titan',
    #     axes=('numclusters', 'ploidy'),
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.parse_titan,
    #     args=(
    #         mgd.InputFile('titan_segs.csv', 'numclusters', 'ploidy', template=segs_template),
    #         mgd.InputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
    #         mgd.InputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
    #         mgd.OutputFile('titan_parsed.csv', 'numclusters', 'ploidy', template=parsed_template),
    #         config['parse_titan'],
    #         sample_id,
    #     ),
    # )
    #
    # workflow.transform(
    #     name='segments_h5',
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.merge_to_h5,
    #     args=(
    #         mgd.InputFile('titan_segs.csv', 'numclusters', 'ploidy', template=segs_template),
    #         mgd.OutputFile(segments),
    #         intervals
    #     ),
    # )
    #
    # workflow.transform(
    #     name='params_h5',
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.merge_to_h5,
    #     args=(
    #         mgd.InputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
    #         mgd.OutputFile(params),
    #         intervals
    #     ),
    # )
    #
    # workflow.transform(
    #     name='markers_h5',
    #     ctx={'mem': global_config['memory']['low'],
    #          'ncpus': 1, 'walltime': '02:00'},
    #     func=tasks.merge_to_h5,
    #     args=(
    #         mgd.InputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
    #         mgd.OutputFile(markers),
    #         intervals
    #     ),
    #     kwargs={'dtype': {'Chr': str}}
    # )

    return workflow
