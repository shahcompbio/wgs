import multiprocessing
import os
import shutil

import pypeliner
import pypeliner.managed as mgd


def run_destruct_local(
        tempdir, tumour_bam, normal_bam,
        sample_id, raw_breakpoints, raw_library,
        reads, destruct_config, refdata_destruct,# logfile,
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

    workflow.subworkflow(
        name='destruct',
        func='destruct.workflow.create_destruct_workflow',
        args=(
            {sample_id: mgd.InputFile(tumour_bam),
             sample_id + 'N': mgd.InputFile(normal_bam)},
            mgd.OutputFile(raw_breakpoints),
            mgd.OutputFile(raw_library),
            mgd.OutputFile(reads),
            destruct_config,
            refdata_destruct
        )
    )

    pyp.run(workflow)

    # pyp_logfile = os.path.join(pipelinedir, 'log', 'latest', 'pipeline.log')
    #
    # if not os.path.exists(logfile):
    #     raise IOError("Couldnt find log file at {}".format(pyp_logfile))
    #
    # shutil.copyfile(pyp_logfile, logfile)
