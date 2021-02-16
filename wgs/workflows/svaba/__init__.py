import pypeliner
import pypeliner.managed as mgd

from wgs.config import config
from wgs.utils import helpers


def create_svaba_workflow(
        tumour_bam,
        normal_bam,
        svaba_vcf,
        reference,
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='run_svaba',
        ctx=helpers.get_default_ctx(
            memory=10,
            walltime='72:00',
            ncpus='8',
            disk=300
        ),
        func='wgs.workflows.svaba.tasks.run_svaba',
        args=(
            mgd.InputFile(tumour_bam),
            mgd.InputFile(normal_bam),
            mgd.TempOutputFile('germline.indel.vcf.gz'),
            mgd.TempOutputFile('germline.sv.vcf.gz'),
            mgd.TempOutputFile('somatic.indel.vcf.gz'),
            mgd.OutputFile(svaba_vcf),
            mgd.TempOutputFile('unfiltered.germline.indel.vcf.gz'),
            mgd.TempOutputFile('unfiltered.germline.sv.vcf.gz'),
            mgd.TempOutputFile('unfiltered.somatic.indel.vcf.gz'),
            mgd.TempOutputFile('unfiltered.somatic.sv.vcf.gz'),
            reference,
            mgd.TempSpace('svaba_tempdir_full')
        ),
        kwargs={
            'ncores': 8,
            'docker_image': config.containers('svaba')
        }
    )

    return workflow
