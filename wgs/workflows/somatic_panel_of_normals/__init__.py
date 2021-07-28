import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


def mutect_workflow(normal_bam, sample_vcf, refdir, single_node):
    chromosomes = config.refdir_data(refdir)['params']['chromosomes']
    paths_refdir = config.refdir_data(refdir)['paths']
    reference = paths_refdir['reference']
    germline_resource = paths_refdir['mutect_panel_of_normals_germline_ref']
    params = config.default_params('variant_calling')

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='generate_intervals',
        func='wgs.workflows.mutect.tasks.generate_intervals',
        ctx=helpers.get_default_ctx(
            memory=5,
            walltime='1:00',
        ),
        ret=mgd.OutputChunks('interval'),
        args=(
            reference,
            chromosomes
        ),
        kwargs={'size': params['split_size']}
    )

    if single_node:
        workflow.transform(
            name='mutect_one_node',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='48:00',
                ncpus=8,
                disk=600
            ),
            axes=('sample_id',),
            func='wgs.workflows.somatic_panel_of_normals.tasks.run_mutect_one_job',
            args=(
                mgd.TempSpace("run_mutect_temp"),
                mgd.OutputFile(sample_vcf),
                reference,
                mgd.InputChunks('interval'),
                mgd.InputFile(normal_bam),
                germline_resource
            ),
        )
    else:
        workflow.transform(
            name='mutect_caller',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='24:00',
            ),
            axes=('interval',),
            func='wgs.workflows.somatic_panel_of_normals.tasks.run_mutect',
            args=(
                mgd.TempOutputFile('mutect.vcf', 'interval'),
                reference,
                mgd.InputInstance('interval'),
                mgd.InputFile(normal_bam),
                germline_resource
            ),
        )

        workflow.transform(
            name='merge_vcfs',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='8:00',
            ),
            func='wgs.workflows.mutect.tasks.merge_vcfs',
            args=(
                mgd.TempInputFile('mutect.vcf', 'interval'),
                mgd.OutputFile(sample_vcf),
                mgd.TempSpace('merge_vcf'),
            ),
        )


def create_somatic_panel_of_normals_workflow(
        samples, normals, pon_vcf, sample_vcfs, refdir, single_node
):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples,
    )

    workflow.subworkflow(
        name="run_mutect_tumour_only",
        func='wgs.workflows.somatic_panel_of_normals.mutect_workflow',
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal.bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.OutputFile('sample.vcf.gz', 'sample_id', fnames=sample_vcfs),
            refdir,
            single_node
        )
    )

    workflow.transform(
        name='create_panel_of_normals',
        ctx=helpers.get_default_ctx(
            memory=15,
            walltime='24:00',
        ),
        func='wgs.workflows.somatic_panel_of_normals.tasks.mutect_panel_of_normals',
        args=(
            mgd.InputFile('sample.vcf.gz', 'sample_id', fnames=sample_vcfs),
            mgd.OutputFile(pon_vcf),
            mgd.TempSpace('pon_create_tempdir')
        ),

    )

    return workflow
