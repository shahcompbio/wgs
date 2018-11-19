import os
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import destruct
from wgs.workflows import lumpy
from wgs.workflows import breakpoint_calling_consensus

def sv_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    destruct_outdir = os.path.join(args['out_dir'], '{sample_id}', 'destruct')
    destruct_breakpoints = os.path.join(destruct_outdir, 'destruct_breakpoints.csv')
    destruct_library = os.path.join(destruct_outdir, 'destruct_library.csv')

    destruct_raw_breakpoints = os.path.join(destruct_outdir, 'destruct_raw_breakpoints.csv')
    destruct_raw_library = os.path.join(destruct_outdir, 'destruct_raw_library.csv')

    destruct_reads = os.path.join(destruct_outdir, 'destruct_reads.csv')
    workflow.subworkflow(
        name='destruct',
        func=destruct.create_destruct_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile(destruct_raw_breakpoints, 'sample_id'),
            mgd.OutputFile(destruct_raw_library, 'sample_id'),
            mgd.OutputFile(destruct_breakpoints, 'sample_id'),
            mgd.OutputFile(destruct_library, 'sample_id'),
            mgd.OutputFile(destruct_reads, 'sample_id'),
            mgd.InputInstance('sample_id'),
            config['globals'],
            config['sv_calling']
        )
    )

    lumpy_outdir = os.path.join(args['out_dir'], '{sample_id}', 'lumpy')
    lumpy_vcf = os.path.join(lumpy_outdir, 'lumpy.vcf')
    workflow.subworkflow(
        name='lumpy',
        func=lumpy.create_lumpy_workflow,
        axes=('sample_id',),
        args=(
            mgd.OutputFile(lumpy_vcf, 'sample_id'),
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

    outdir = os.path.join(args['out_dir'], '{sample_id}')
    parsed_csv = os.path.join(outdir, 'filtered_consensus_calls.csv')
    workflow.subworkflow(
        name="consensus_calling",
        func=breakpoint_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile(destruct_breakpoints, 'sample_id'),
            mgd.InputFile(lumpy_vcf, 'sample_id'),
            mgd.OutputFile(parsed_csv, 'sample_id'),
            config['globals'],
            config['sv_calling'],
            mgd.InputInstance('sample_id')
        ),
    )

    pyp.run(workflow)
