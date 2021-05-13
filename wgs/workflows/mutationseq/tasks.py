'''
Created on Feb 21, 2018

@author: dgrewal
'''
import pypeliner
import pysam

from wgs.workflows.mutationseq.scripts import PlotSingleSample


def generate_intervals(ref, chromosomes, size=1000000):
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    names = fasta.references

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size) + 1)
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    return intervals


def run_museqportrait(
        infile, out_pdf, out_txt, museqportrait_log,
        single_mode, thousand_genomes=None, dbsnp=None, germline_refdata=None, germline_plot_threshold=0.5
):
    """
    Run museqportrait script on the input VCF file

    :param infile:
    :type infile: str
    :param out_pdf:
    :type out_pdf: str
    :param out_txt:
    :type out_txt: str
    :param museqportrait_log:
    :type museqportrait_log: str
    :param single_mode:
    :type single_mode: bool
    :param plot_params:
    :type plot_params: dict
    :param databases:
    :type databases: dict
    :return:
    :rtype:
    """

    if single_mode:
        plt_ss = PlotSingleSample(
            infile,
            thousand_genomes,
            dbsnp,
            germline_refdata,
            out_pdf,
            germline_plot_threshold
        )

        plt_ss.main()

        # touch the txt file to avoid pypeliner errors
        open(out_txt, 'w').close()
        open(museqportrait_log, 'w').close()

    else:
        cmd = ['museqportrait', '--log', museqportrait_log, '--output-pdf',
               out_pdf, '--output-txt', out_txt, infile]
        pypeliner.commandline.execute(*cmd)


def fix_museq_vcf(infile, output):
    with open(infile, 'rt') as reader, open(output, 'wt') as writer:
        for line in reader:
            if line.startswith('#'):
                line = line.replace('##FORMAT=<ID=PL,Number=3', '##FORMAT=<ID=PL,Number=G')
            writer.write(line)
