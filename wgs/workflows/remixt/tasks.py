import logging
import multiprocessing
import os

import pandas as pd
import pypeliner
import pypeliner.managed as mgd


def filter_destruct_breakpoints(breakpoints, filtered_breakpoints, min_num_reads):
    breakpoints = pd.read_csv(breakpoints, sep='\t')

    breakpoints = breakpoints[breakpoints['num_reads'] >= min_num_reads]

    breakpoints = breakpoints[[
        'prediction_id',
        'chromosome_1',
        'strand_1',
        'position_1',
        'chromosome_2',
        'strand_2',
        'position_2',
    ]]

    breakpoints.to_csv(filtered_breakpoints, sep='\t', index=False)


def run_remixt_local(
        tempdir, breakpoints, tumour_bam, normal_bam,
        sample_id, remixt_results_filename, remixt_raw_dir,
        remixt_config, remixt_refdata,
        ncpus=None
):
    pipelinedir = os.path.join(tempdir, 'pipeline')
    tmpdir = os.path.join(tempdir, 'tmp')

    if not ncpus:
        ncpus = multiprocessing.cpu_count()

    config = {'pipelinedir': pipelinedir, 'tmpdir': tmpdir,
              'submit': 'local', 'maxjobs': ncpus,
              'loglevel': 'DEBUG'}

    pyp = pypeliner.app.Pypeline(config=config)
    workflow = pypeliner.workflow.Workflow()

    logging.getLogger().setLevel(logging.DEBUG)

    workflow.subworkflow(
        name='remixt',
        func="remixt.workflow.create_remixt_bam_workflow",
        args=(
            mgd.InputFile(breakpoints),
            {sample_id: mgd.InputFile(tumour_bam),
             sample_id + 'N': mgd.InputFile(normal_bam)},
            {sample_id: mgd.OutputFile(remixt_results_filename)},
            remixt_raw_dir,
            remixt_config,
            remixt_refdata,
        ),
        kwargs={
            'normal_id': sample_id + 'N',
        }
    )

    pyp.run(workflow)
