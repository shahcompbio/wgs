import pypeliner
import pypeliner.managed as mgd

import tasks


def create_ichorcna_workflow(
        tumour_bam, normal_panel, segments,
        params, depth, config, global_config,
        sample_id, plots_tar
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='hmmcopy_readcounter',
        func=tasks.hmmcopy_readcounter,
        ctx={
            'mem': global_config['memory']['low'],
            'ncpus': 1,
        },
        args=(
            mgd.InputFile(tumour_bam),
            mgd.TempOutputFile('infile.wig'),
            config,
        )
    )

    workflow.transform(
        name='run_ichorcna',
        func=tasks.run_ichorcna,
        ctx={
            'mem': global_config['memory']['med'],
            'ncpus': 1,
        },
        args=(
            mgd.TempInputFile('infile.wig'),
            mgd.InputFile(normal_panel),
            mgd.OutputFile(segments),
            mgd.OutputFile(params),
            mgd.OutputFile(depth),
            config['centromere'],
            config['gc_wig'],
            config['map_wig'],
            sample_id,
            mgd.TempSpace('plots_dir'),
            mgd.OutputFile(plots_tar)
        ),
        kwargs=dict(txnE=config['txnE'], chromosomes=config['chromosomes'])
    )

    return workflow
