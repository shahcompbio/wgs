import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import alignment


def realign_bam_files(
        input, output, metrics_output,
        metrics_tar, refdir,
        single_node=False,
        ignore_bamtofastq_exception=False,
        picard_mem=8
):
    outputs_tdf = output + '.tdf'

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='bam_to_fastq',
        ctx=helpers.get_default_ctx(
            walltime='96:00',
            disk=500
        ),
        func="wgs.workflows.realignment.tasks.split_by_rg",
        args=(
            mgd.InputFile(input),
            mgd.TempOutputFile("inputdata_read1.fastq.gz",  "readgroup"),
            mgd.TempOutputFile("inputdata_read2.fastq.gz",  "readgroup", axes_origin=[]),
            mgd.TempSpace("bamtofastq"),
            ignore_bamtofastq_exception
        )
    )

    workflow.transform(
        name='get_sample_info',
        func="wgs.workflows.realignment.tasks.get_read_group",
        ret=mgd.TempOutputObj('sample_info'),
        args=(
            mgd.InputFile(input),
        )
    )

    workflow.subworkflow(
        name='align_samples',
        func=alignment.align_samples,
        args=(
            mgd.TempInputFile("inputdata_read1.fastq.gz", "readgroup"),
            mgd.TempInputFile("inputdata_read2.fastq.gz", "readgroup"),
            mgd.OutputFile(output, extensions=['.bai']),
            mgd.OutputFile(metrics_output, extensions=['.yaml']),
            mgd.OutputFile(metrics_tar),
            mgd.OutputFile(outputs_tdf),
            mgd.TempInputObj('sample_info'),
            refdir
        ),
        kwargs={
            'single_node': single_node,
            'picard_mem': picard_mem
        }
    )

    return workflow
