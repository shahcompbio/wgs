import os
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import mutationseq
from wgs.workflows import strelka
from wgs.workflows import vcf_annotation
from wgs.workflows import variant_calling_consensus


def variant_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])
    output_dir = os.path.join(args['out_dir'], '{sample_id}')

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}

    museq_vcf = os.path.join(output_dir, '{sample_id}', 'museq_paired_annotated.vcf')
    museq_ss_vcf = os.path.join(output_dir, '{sample_id}', 'museq_single_annotated.vcf')
    strelka_snv_vcf = os.path.join(output_dir, '{sample_id}', 'strelka_snv_annotated.vcf')
    strelka_indel_vcf = os.path.join(output_dir, '{sample_id}', 'strelka_indel_annotated.vcf')

    parsed_csv = os.path.join(output_dir, '{sample_id}', 'allcalls.csv')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    museqportrait_pdf = os.path.join(output_dir, 'paired_museqportrait.pdf')
    museqportrait_txt = os.path.join(output_dir, 'paired_museqportrait.txt')
    workflow.subworkflow(
        name="mutationseq_paired",
        func=mutationseq.create_museq_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile(museqportrait_pdf, 'sample_id'),
            mgd.OutputFile(museqportrait_txt, 'sample_id'),
            config['globals'],
            config['variant_calling'],
            mgd.InputInstance('sample_id')
        ),
        kwargs={
            'tumour_bam': mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                                        extensions=['.bai'], axes_origin=[]),
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
        }
    )

    museqportrait_pdf = os.path.join(output_dir, 'single_museqportrait.pdf')
    museqportrait_txt = os.path.join(output_dir, 'single_museqportrait.txt')
    workflow.subworkflow(
        name="mutationseq_single",
        func=mutationseq.create_museq_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempOutputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile(museqportrait_pdf, 'sample_id'),
            mgd.OutputFile(museqportrait_txt, 'sample_id'),
            config['globals'],
            config['variant_calling'],
            mgd.InputInstance('sample_id')
        ),
        kwargs={
            'tumour_bam': None,
            'normal_bam': mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                                        extensions=['.bai'], axes_origin=[]),
        }
    )

    workflow.subworkflow(
        name="strelka",
        func=strelka.create_strelka_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            config['reference'],
            mgd.TempOutputFile('strelka_indel.vcf.gz', 'sample_id'),
            mgd.TempOutputFile('strelka_snv.vcf.gz', 'sample_id'),
            config['globals'],
            config['variant_calling'],
        ),
    )

    workflow.subworkflow(
        name="annotate_paired_museq",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_snv_ann.vcf.gz', 'sample_id',
                               extensions=['.csi','.tbi'], template=museq_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
    )

    workflow.subworkflow(
        name="annotate_germline_museq",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("museq_germlines.vcf.gz", 'sample_id'),
            mgd.OutputFile('museq_germlines_ann.vcf.gz', 'sample_id',
                               extensions=['.csi','.tbi'], template=museq_ss_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
    )

    workflow.subworkflow(
        name="annotate_strelka",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_snv.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_snv_ann.vcf.gz', 'sample_id',
                               extensions=['.csi','.tbi'], template=strelka_snv_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
    )

    workflow.subworkflow(
        name="annotate_strelka_indel",
        func=vcf_annotation.create_annotation_workflow,
        axes=('sample_id',),
        args=(
            mgd.TempInputFile("strelka_indel.vcf.gz", 'sample_id'),
            mgd.OutputFile('strelka_indel_ann.vcf.gz', 'sample_id',
                               extensions=['.csi','.tbi'], template=strelka_indel_vcf),
            config['globals'],
            config['variant_calling']['annotation_params'],
        ),
    )

    workflow.subworkflow(
        name="consensus_calling",
        func=variant_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("museq_germlines_ann.vcf.gz", 'sample_id', template=museq_ss_vcf),
            mgd.InputFile("museq_snv_ann.vcf.gz", 'sample_id', template=museq_vcf),
            mgd.InputFile("strelka_snv_ann.vcf.gz", 'sample_id', template=strelka_snv_vcf),
            mgd.InputFile("strelka_indel_ann.vcf.gz", 'sample_id', template=strelka_indel_vcf),
            mgd.OutputFile(parsed_csv, 'sample_id'),
            config['globals'],
            config['variant_calling'],
        ),
    )

    pyp.run(workflow)
