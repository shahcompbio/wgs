import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers

from wgs.config import config

def create_hmmcopy_workflow(
        bam_file, sample_id, bias_pdf, correction_pdf,
        hmmcopy_pdf, hmmcopy_table, pygenes_table,
        chromosomes, map_wig, gc_wig, pygenes_gtf,
):


    cn_params = config.default_params()['copynumber_calling']

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='hmmcopy_readcounter',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='2:00', ),
        func='wgs.workflows.hmmcopy.tasks.hmmcopy_readcounter',
        args=(
            mgd.InputFile(bam_file, extensions=['.bai']),
            mgd.TempOutputFile('infile.wig'),
            chromosomes,
            cn_params['readcounter'],
        )
    )

    workflow.transform(
        name='calc_corr',
        func='wgs.workflows.hmmcopy.tasks.calc_corr',
        args=(
            mgd.TempInputFile('infile.wig'),
            mgd.TempOutputFile('infile_copy.txt'),
            mgd.TempOutputFile('infile_copy.obj'),
            gc_wig,
            map_wig,
            cn_params['map_cutoff'],
        ),
        kwargs={'docker_image': config.containers('hmmcopy')}
    )

    workflow.transform(
        name='run_hmmcopy',
        func='wgs.workflows.hmmcopy.tasks.run_hmmcopy',
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('infile_copy.txt'),
            mgd.TempOutputFile('hmmcopy_res.obj'),
            mgd.TempOutputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(hmmcopy_table),
            sample_id,
            cn_params['hmmcopy_params'],
        ),
        kwargs={'docker_image': config.containers('hmmcopy')}
    )

    workflow.transform(
        name='plot_hmm',
        func='wgs.workflows.hmmcopy.tasks.plot_hmm',
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('hmmcopy_res.obj'),
            mgd.TempSpace('correction_plots_dir'),
            mgd.TempSpace('hmmcopy_plots_dir'),
            mgd.OutputFile(bias_pdf),
            mgd.OutputFile(correction_pdf),
            mgd.OutputFile(hmmcopy_pdf),
        ),
        kwargs={'docker_image': config.containers('hmmcopy')}
    )

    workflow.transform(
        name='annot_hmm',
        func='wgs.workflows.hmmcopy.tasks.annot_hmm',
        args=(
            mgd.TempInputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(pygenes_table),
            pygenes_gtf,
        )
    )

    return workflow
