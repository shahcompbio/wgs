import logging
import multiprocessing
import os

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def run_destruct_local(
        tempdir, tumour_bam, normal_bam,
        sample_id, raw_breakpoints, raw_library,
        reads, destruct_config, refdata_destruct,
        ncpus=None, docker_image=None
):
    pipelinedir = os.path.join(tempdir, 'pipeline')
    tmpdir = os.path.join(tempdir, 'tmp')

    if not ncpus:
        ncpus = multiprocessing.cpu_count()

    config = {'pipelinedir': pipelinedir, 'tmpdir': tmpdir,
              'submit': 'local', 'maxjobs': ncpus,
              'loglevel': 'DEBUG'}

    pyp = pypeliner.app.Pypeline(config=config)

    context_config = pypeliner.helpers.GlobalState.get_all()['context_config']
    pypeliner.helpers.GlobalState.set('context_config', context_config)

    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': docker_image})

    logging.getLogger().setLevel(logging.DEBUG)

    workflow.subworkflow(
        name='destruct_local_in_job',
        func='destruct.workflow.create_destruct_workflow',
        ctx={'docker_image': docker_image},
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


def reheader_reads(infile, outfile):
    with helpers.GetFileHandle(infile, 'rt') as indata, helpers.GetFileHandle(outfile, 'wt') as outdata:
        line_one = indata.readline()

        header = line_one.split('\t')

        assert len(header) == 8

        assert not header[0] == 'prediction_id'

        header = ['prediction_id', 'library_id', 'fragment_id', 'read_end', 'seq', 'qual', 'comment', 'filtered']

        header = '\t'.join(header) + '\n'

        outdata.write(header)
        outdata.write(line_one)

        for line in indata:
            outdata.write(line)
