import os

import pypeliner
import pypeliner.managed as mgd
from  wgs.workflows.hmmcopy import tasks
from wgs.utils import helpers


def create_hmmcopy_workflow(
        bam_file, out_dir, global_config, config,
        sample_id, bias_pdf, correction_pdf, hmmcopy_pdf,
        hmmcopy_table, pygenes_table
):

    workflow = pypeliner.workflow.Workflow()


    workflow.transform(
        name='hmmcopy_readcounter',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['low'],
            walltime='2:00', ),
        func=tasks.hmmcopy_readcounter,
        args=(
            mgd.InputFile(bam_file, extensions=['.bai']),
            mgd.TempOutputFile('infile.wig'),
            config,
        )
    )

    workflow.transform(
        name='calc_corr',
        func=tasks.calc_corr,
        args=(
            mgd.TempInputFile('infile.wig'),
            mgd.TempOutputFile('infile_copy.txt'),
            mgd.TempOutputFile('infile_copy.obj'),
            config,
        ),
        kwargs={'docker_image': config['docker']['hmmcopy']}
    )

    workflow.transform(
        name='run_hmmcopy',
        func=tasks.run_hmmcopy,
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('infile_copy.txt'),
            mgd.TempOutputFile('hmmcopy_res.obj'),
            mgd.TempOutputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(tumour_table_out),
            sample_id,
            config,
        ),
        kwargs={'docker_image': config['docker']['hmmcopy']}
    )

    workflow.transform(
        name='plot_hmm',
        func=tasks.plot_hmm,
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('hmmcopy_res.obj'),
            mgd.TempSpace('correction_plots_dir'),
            mgd.TempSpace('hmmcopy_plots_dir'),
            mgd.OutputFile(bias_plots_pdf),
            mgd.OutputFile(correction_plots_pdf),
            mgd.OutputFile(hmmcopy_plots_pdf),
        ),
        kwargs={'docker_image': config['docker']['hmmcopy']}
    )

    workflow.transform(
        name='annot_hmm',
        func=tasks.annot_hmm,
        args=(
            mgd.TempInputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(pygene_outfile),
            config,
        )
    )

    return workflow
