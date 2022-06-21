'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.config import config
from wgs.utils import helpers


class MissingInputs(Exception):
    pass


def create_remixt_workflow(
        tumour_path,
        normal_path,
        breakpoints,
        sample_id,
        remixt_results_filename,
        remixt_brk_cn_csv,
        remixt_cn_csv,
        remixt_minor_modes_csv,
        remixt_mix_csv,
        remixt_read_depth_csv,
        remixt_stats_csv,
        remixt_refdata,
        remixt_raw_dir,
        reference,
        chromosomes,
        single_node=False,
):
    params = config.default_params('copynumber_calling')['remixt']

    workflow = pypeliner.workflow.Workflow()

    if chromosomes[0].startswith('chr'):
        chr_prefix = 'chr'
        chromosomes = [v.replace('chr', '') for v in chromosomes]
    else:
        chr_prefix = ''

    remixt_config = {
        'genome_fasta_template': reference,
        'genome_fai_template': reference + '.fai',
        'max_copy_number': 10,
        'max_ploidy': 7,
        'chromosomes': chromosomes,
        'chr_prefix': chr_prefix
    }

    if breakpoints is None:
        raise MissingInputs('Remixt requires breakpoints')

    workflow.transform(
        name='filter_breakpoints',
        func='wgs.workflows.remixt.tasks.filter_destruct_breakpoints',
        ctx=helpers.get_default_ctx(
            memory=4,
            walltime='4:00'),
        args=(
            mgd.InputFile(breakpoints),
            mgd.TempOutputFile('filtered_breakpoints.csv'),
            params['min_num_reads']
        )
    )

    if single_node:
        workflow.transform(
            name='remixt',
            func='wgs.workflows.remixt.tasks.run_remixt_local',
            ctx=helpers.get_default_ctx(
                memory=15,
                walltime='120:00',
                ncpus=8),
            args=(
                mgd.TempSpace("remixt_temp"),
                mgd.TempInputFile('filtered_breakpoints.csv'),
                mgd.InputFile(tumour_path, extensions=['.bai']),
                mgd.InputFile(normal_path, extensions=['.bai']),
                sample_id,
                mgd.OutputFile(remixt_results_filename),
                remixt_raw_dir,
                remixt_config,
                remixt_refdata,
            ),
        )
    else:
        workflow.subworkflow(
            name='remixt',
            func="remixt.workflow.create_remixt_bam_workflow",
            ctx={'walltime': '48:00'},
            args=(
                mgd.TempInputFile('filtered_breakpoints.csv'),
                {sample_id: mgd.InputFile(tumour_path, extensions=['.bai']),
                 sample_id + 'N': mgd.InputFile(normal_path, extensions=['.bai'])},
                {sample_id: mgd.OutputFile(remixt_results_filename)},
                remixt_raw_dir,
                remixt_config,
                remixt_refdata,
            ),
            kwargs={
                'normal_id': sample_id + 'N',
            }
        )

    workflow.transform(
        name='parse_remixt',
        func='wgs.workflows.remixt.tasks.parse_remixt_file',
        args=(
            mgd.InputFile(remixt_results_filename),
            [
                mgd.OutputFile(remixt_brk_cn_csv, extensions=['.yaml']),
                mgd.OutputFile(remixt_cn_csv, extensions=['.yaml']),
                mgd.OutputFile(remixt_minor_modes_csv, extensions=['.yaml']),
                mgd.OutputFile(remixt_mix_csv, extensions=['.yaml']),
                mgd.OutputFile(remixt_read_depth_csv, extensions=['.yaml']),
                mgd.OutputFile(remixt_stats_csv, extensions=['.yaml']),
            ],
            [
                '/brk_cn', '/cn', '/minor_modes', '/mix', '/read_depth', '/stats'
            ],
            mgd.TempSpace('tempdir_parse')
        )
    )

    return workflow
