import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import alignment


def realign_bam_file(input, output, outdir, config, single_node=False):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='bam_to_fastq',
        ctx=helpers.get_default_ctx(
            walltime='72:00',
            disk=500
        ),
        func="wgs.workflows.realignment.tasks.split_by_rg",
        args=(
            mgd.InputFile(input),
            mgd.TempOutputFile("inputdata_read1.fastq", "readgroup"),
            mgd.TempOutputFile("inputdata_read2.fastq", "readgroup", axes_origin=[]),
            mgd.TempSpace("bamtofastq")
        )
    )
    workflow.subworkflow(
        name='align_samples',
        func=alignment.align_sample,
        axes=('readgroup',),
        args=(
            config,
            mgd.TempInputFile("inputdata_read1.fastq", "readgroup"),
            mgd.TempInputFile("inputdata_read2.fastq", "readgroup"),
            mgd.TempOutputFile('aligned_lanes.bam', 'readgroup'),
            outdir,
            mgd.InputInstance("readgroup"),
        ),
        kwargs={'single_node': single_node}
    )

    workflow.transform(
        name='merge_tumour_lanes',
        ctx=helpers.get_default_ctx(
            walltime='96:00',
            disk=500
        ),
        func="wgs.workflows.alignment.tasks.merge_bams",
        args=(
            mgd.TempInputFile('aligned_lanes.bam', 'readgroup'),
            mgd.OutputFile(output, extensions=['.bai']),
        ),
        kwargs={
            'picard_docker_image': config['docker']['picard'],
            'samtools_docker_image': config['docker']['samtools']
        }
    )

    return workflow
