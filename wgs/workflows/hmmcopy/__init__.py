import os

import pypeliner
import pypeliner.managed as mgd
from  wgs.workflows.hmmcopy import tasks
from wgs.utils import helpers


def create_hmmcopy_workflow(bam_file, out_dir, global_config, config, sample_id):
    bias_plots_pdf = os.path.join(out_dir, 'plots', '{}_bias.pdf'.format(sample_id))
    correction_plots_pdf = os.path.join(out_dir, 'plots', '{}_correction.pdf'.format(sample_id))
    hmmcopy_plots_pdf = os.path.join(out_dir, 'plots', '{}_hmmcopy.pdf'.format(sample_id))
    tumour_table_out = os.path.join(out_dir, '{}_tumour_correctreads_with_state.txt'.format(sample_id))
    pygene_outfile = os.path.join(out_dir, '{}_hmmcopy.seg.pygenes'.format(sample_id))

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
