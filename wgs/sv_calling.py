import os
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import destruct
from wgs.workflows import lumpy
from wgs.workflows import breakpoint_calling_consensus


def call_breakpoints(
        samples, config, tumours, normals, destruct_raw_breakpoints,
        destruct_raw_library, destruct_breakpoints, destruct_library,
        destruct_reads, lumpy_vcf, consensus_calls
):
    destruct_raw_breakpoints = dict([(sampid, destruct_raw_breakpoints[sampid])
                                     for sampid in samples])
    destruct_raw_library = dict([(sampid, destruct_raw_library[sampid])
                                for sampid in samples])
    destruct_breakpoints = dict([(sampid, destruct_breakpoints[sampid])
                                for sampid in samples])
    destruct_library = dict([(sampid, destruct_library[sampid])
                            for sampid in samples])
    destruct_reads = dict([(sampid, destruct_reads[sampid])
                          for sampid in samples])
    lumpy_vcf = dict([(sampid, lumpy_vcf[sampid])
                     for sampid in samples])
    consensus_calls = dict([(sampid, consensus_calls[sampid])
                           for sampid in samples])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='destruct',
        func=destruct.create_destruct_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('destruct_raw_breakpoints', 'sample_id', fnames=destruct_raw_breakpoints),
            mgd.OutputFile('destruct_raw_library', 'sample_id', fnames=destruct_raw_library),
            mgd.OutputFile('destruct_breakpoints', 'sample_id', fnames=destruct_breakpoints),
            mgd.OutputFile('destruct_library', 'sample_id', fnames=destruct_library),
            mgd.OutputFile('destruct_reads', 'sample_id', fnames=destruct_reads),
            mgd.InputInstance('sample_id'),
            config['globals'],
            config['sv_calling']
        )
    )

    workflow.subworkflow(
        name='lumpy',
        func=lumpy.create_lumpy_workflow,
        axes=('sample_id',),
        args=(
            mgd.OutputFile('lumpy_vcf', 'sample_id', fnames=lumpy_vcf),
            config['globals'],
            config['sv_calling']
        ),
        kwargs={
            'tumour_bam': mgd.InputFile(
                "tumour.bam", 'sample_id', fnames=tumours,
                extensions=['.bai'], axes_origin=[]),
            'normal_bam': mgd.InputFile(
                "normal.bam", 'sample_id', fnames=normals,
                extensions=['.bai'], axes_origin=[]),
        }
    )

    workflow.subworkflow(
        name="consensus_calling",
        func=breakpoint_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('destruct_breakpoints', 'sample_id', fnames=destruct_breakpoints),
            mgd.InputFile('lumpy_vcf', 'sample_id', fnames=lumpy_vcf),
            mgd.OutputFile('consensus_calls', 'sample_id', fnames=consensus_calls),
            config['globals'],
            config['sv_calling'],
            mgd.InputInstance('sample_id')
        ),
    )

    return workflow


def sv_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}

    sv_outdir = os.path.join(args['out_dir'], 'breakpoints', '{sample_id}')
    destruct_breakpoints = os.path.join(sv_outdir, 'destruct_breakpoints.csv')
    destruct_library = os.path.join(sv_outdir, 'destruct_library.csv')
    destruct_raw_breakpoints = os.path.join(sv_outdir, 'destruct_raw_breakpoints.csv')
    destruct_raw_library = os.path.join(sv_outdir, 'destruct_raw_library.csv')
    destruct_reads = os.path.join(sv_outdir, 'destruct_reads.csv')
    lumpy_vcf = os.path.join(sv_outdir, 'lumpy.vcf')
    parsed_csv = os.path.join(sv_outdir, 'filtered_consensus_calls.csv')

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="call_breakpoints",
        func=call_breakpoints,
        args=(
            samples,
            config,
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('destruct_raw_breakpoints', 'sample_id', template=destruct_raw_breakpoints, axes_origin=[]),
            mgd.OutputFile('destruct_raw_library', 'sample_id', template=destruct_raw_library, axes_origin=[]),
            mgd.OutputFile('destruct_breakpoints', 'sample_id', template=destruct_breakpoints, axes_origin=[]),
            mgd.OutputFile('destruct_library', 'sample_id', template=destruct_library, axes_origin=[]),
            mgd.OutputFile('destruct_reads', 'sample_id', template=destruct_reads, axes_origin=[]),
            mgd.OutputFile('lumpy_vcf', 'sample_id', template=lumpy_vcf, axes_origin=[]),
            mgd.OutputFile('parsed_csv', 'sample_id', template=parsed_csv, axes_origin=[])
        )
    )

    pyp.run(workflow)
