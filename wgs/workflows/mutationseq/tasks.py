'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pysam

from scripts import PlotSingleSample


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size)+1)
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def run_museqportrait(infile, out_pdf, out_txt, museqportrait_log, single_mode, config,
                      docker_image=None):
    '''
    Run museqportrait script on the input VCF file

    :param infile: temporary input VCF file
    :param out_dir: temporary output VCF file
    :param museqportrait_log: path to the log file
    '''

    if single_mode:
        plt_ss = PlotSingleSample(
            infile,
            config["thousandgen_params"]["db"],
            config["dbsnp_params"]["db"],
            config['refdata_single_sample'],
            out_pdf,
            config['threshold']
        )

        plt_ss.main()

        # touch the txt file to avoid pypeliner errors
        open(out_txt, 'w').close()
        open(museqportrait_log, 'w').close()

    else:
        cmd = ['museqportrait', '--log', museqportrait_log, '--output-pdf',
               out_pdf, '--output-txt', out_txt, infile]
        pypeliner.commandline.execute(*cmd, docker_image=docker_image)
