import os

import pypeliner
from wgs.utils import helpers
from wgs.utils import pdfutils
from wgs.workflows.hmmcopy.scripts import PygeneAnnotation
from wgs.workflows.hmmcopy.scripts import ReadCounter

scripts = os.path.join(
    os.path.realpath(os.path.dirname(__file__)),
    'scripts'
)


def hmmcopy_readcounter(input_bam, output_wig, chromosomes, config):
    rc = ReadCounter(
        input_bam, output_wig, config['w'],
        chromosomes, config['q'],
    )
    rc.main()


def calc_corr(input_wig, output_file, output_obj, gc_wig, map_wig, map_cutoff):
    cmd = [
        'hmmcopy_correct_reads.R',
        input_wig,
        gc_wig,
        map_wig,
        map_cutoff,
        output_file,
        output_obj,
    ]

    pypeliner.commandline.execute(*cmd)


def run_hmmcopy(
        tumour_copy,
        tumour_table,
        output_obj,
        output_segments,
        tumour_table_out,
        sample_id,
        hmmcopy_params,
):
    args = {key: 'NULL' if value is None else value
            for key, value in hmmcopy_params.items()}

    param_list = (
        'm',
        'mu',
        'kappa',
        'e',
        'S',
        'strength',
        'lambda',
        'nu',
        'eta',
        'gamma',
    )

    params = [args[param] for param in param_list]

    cmd = [
              'hmmcopy.R',
              tumour_copy,
              tumour_table,
              args['normal_copy'],
              args['normal_table'],
              output_segments,
              output_obj,
              sample_id,
              args['normal_table_out'],
              tumour_table_out,
          ] + params

    pypeliner.commandline.execute(*cmd)


def plot_hmm(
        tumour_copy,
        hmmcopy_res,
        correction_plots_dir,
        hmmcopy_plots_dir,
        bias_pdf,
        correction_pdf,
        hmmcopy_pdf,
        chromosomes
):
    helpers.makedirs(correction_plots_dir)
    helpers.makedirs(hmmcopy_plots_dir)

    cmd = [
        'plot_hmmcopy.R',
        tumour_copy,
        hmmcopy_res,
        correction_plots_dir,
        bias_pdf,
        hmmcopy_plots_dir,
    ]

    pypeliner.commandline.execute(*cmd)

    correction_pdfs = [os.path.join(correction_plots_dir, f)
                       for f in os.listdir(correction_plots_dir) if f.endswith('.jpg')]
    pdfutils.merge_jpgs(correction_pdfs, correction_pdf)

    all_hmmcopy_pdfs = [os.path.join(hmmcopy_plots_dir, pdf)
                        for pdf in os.listdir(hmmcopy_plots_dir)]
    # just some sorting
    pdfs = [os.path.join(hmmcopy_plots_dir, 'chr_{}.pdf'.format(chrom))
                  for chrom in chromosomes]
    all_hmmcopy_pdfs = [v for v in pdfs if v in all_hmmcopy_pdfs]
    all_hmmcopy_pdfs += list(set(all_hmmcopy_pdfs) - set(pdfs))

    pdfutils.merge_pdfs(pdfs, hmmcopy_pdf)


def annot_hmm(input_segments, output_file, pygenes_gtf):
    annotator = PygeneAnnotation(input_segments, output_file, gtf=pygenes_gtf)
    annotator.write_output()
