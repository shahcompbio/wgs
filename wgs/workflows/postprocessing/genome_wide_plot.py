import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf

from wgs_qc_utils.plotter import variant_plotting, coverage_plotting, \
    ideogram_plotting, remixt_plotting, roh_plotting, snv_cn, \
    titan_plotting, gene_annotation_plotting
from wgs_qc_utils.reader import read_remixt, read_roh, read_variant_calls, \
    read_titan, read_coverage, parse_snv_cn
from wgs_qc_utils.reader.ideogram import read_ideogram
import logging


def plot_chrom_on_axes(remixt, titan, roh, germline_calls, somatic_calls,
                       tumour_coverage, normal_coverage, breakpoints, vaf_data,
                       ideogram, chrom, axes):
    """
    plot one chromosome onto a set of axes
    :param copy_number: read-in copy number data
    :param roh: read-in roh data
    :param germline_calls: read-in germline cals
    :param somatic_calls: read-in somatic calls
    :param tumour_coverage: read-in tumour coverage
    :param normal_coverage: read-in normal coverage
    :param chrom: chromosome to plot
    :param axes: set of axes to plot on
    :param remixt: T/F: true if copy number data is Remixt, false if Titan
    :return: axes with plots added
    """
    prepped_remxit = read_remixt.prepare_at_chrom(remixt, chrom)
    prepped_titan = read_titan.prepare_at_chrom(titan, chrom)
    prepped_roh = read_roh.prepare_at_chrom(roh, chrom)
    prepped_somatic_calls = read_variant_calls.prepare_at_chrom(somatic_calls, chrom, n_bins=2000)
    prepped_germline_calls = read_variant_calls.prepare_at_chrom(germline_calls, chrom, n_bins=2000)
    prepped_tumour_coverage = read_coverage.prepare_at_chrom(tumour_coverage, chrom)
    prepped_normal_coverage = read_coverage.prepare_at_chrom(normal_coverage, chrom)
    prepped_breakpoints = read_variant_calls.prepare_at_chrom(breakpoints, chrom, n_bins=2000)
    prepped_snv_cn = parse_snv_cn.prepare_at_chrom(vaf_data, chrom)
    prepped_ideogram = read_ideogram.prepare_at_chrom(ideogram, chrom)

    # coverage_ylim_max = prepped_tumour_coverage.coverage.max() + 10
    # coverage_ylim_min = prepped_normal_coverage.coverage.min() - 10
    coverage_ylim_min = 0
    coverage_ylim_max = 150
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_max = 250
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_min = 0
    anno_genes = gene_annotation_plotting.get_gene_annotation_data(chrom)

    chrom_max = prepped_ideogram.start.max()
    axes[0] = variant_plotting.plot_fill(prepped_somatic_calls.location, prepped_somatic_calls.n_events,
                                         axes[0], "Somatic \n Alteration \n Frequency", chrom_max)
    axes[1] = variant_plotting.plot_fill(prepped_germline_calls.location, prepped_germline_calls.n_events,
                                         axes[1], "Germline \n Alteration \n Frequency", chrom_max)
    axes[2] = coverage_plotting.plot(prepped_tumour_coverage.start, prepped_tumour_coverage.coverage,
                                     coverage_ylim_min, coverage_ylim_max,
                                     axes[2], "Tumour \n Coverage", chrom_max)
    axes[3] = coverage_plotting.plot(prepped_normal_coverage.start, prepped_normal_coverage.coverage,
                                     coverage_ylim_min, coverage_ylim_max,
                                     axes[3], "Normal \n Coverage", chrom_max)

    axes[4] = variant_plotting.plot_bar(prepped_breakpoints.location, prepped_breakpoints.n_events,
                                        axes[4], "Breakpoint \n Frequency", chrom_max)

    axes[5] = remixt_plotting.plot(prepped_remxit.start, prepped_remxit.major_raw,
                                    prepped_remxit.minor_raw, axes[5], chrom_max)

    axes[5] = snv_cn.plot_scatter(prepped_snv_cn.pos,prepped_snv_cn.frac_cn, axes[5])

    axes[6] = titan_plotting.plot(prepped_titan.Position,prepped_titan.LogRatio,
                                  prepped_titan.color, axes[6], chrom_max, anno_genes=anno_genes)

    axes[7] = roh_plotting.plot(prepped_roh.pos, prepped_roh.state, axes[7], chrom_max)

    axes[8] = ideogram_plotting.plot(prepped_ideogram, axes[8])

    axes[14] = snv_cn.plot_hist(prepped_snv_cn.frac_cn, axes[14])

    axes[23] = remixt_plotting.add_remixt_legend(axes[23])
    axes[15] = titan_plotting.add_titan_legend(axes[15])

    if not anno_genes.empty:
        axes[24] = gene_annotation_plotting.add_gene_annotation_legend(anno_genes, axes[24])

    return axes


def rasturize_axes(axes):
    """
    rasterize axes
    :param axes: list of axes to rasterize
    :return: rasterized axes
    """
    for ax in axes:
        ax.set_rasterized(True)
    return axes


def genome_wide_plot(remixt, remixt_label, titan, roh, germline_calls, somatic_calls,
                     tumour_coverage, normal_coverage, breakpoints, chromosomes, pdf):

    """
    make a genome wide plot
    :param copy_number: copy number data
    :param roh: roh data
    :param germline_calls: germline data
    :param somatic_calls: somatic data
    :param tumour_coverage: tumour coverage data
    :param normal_coverage: normal coverage data
    :param sample_label: label for sample plotted
    :param plot_name: name for plot pdf
    :param remixt: remixt data
    """

    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf)

    remixt = read_remixt.read(remixt, remixt_label)

    titan = read_titan.read(titan)

    roh = read_roh.read(roh)

    germline_calls = read_variant_calls.read(germline_calls)
    somatic_calls = read_variant_calls.read_consensus_csv(somatic_calls)

    tumour_coverage = read_coverage.read(tumour_coverage)
    normal_coverage = read_coverage.read(normal_coverage)

    breakpoints = read_variant_calls.read_svs(breakpoints)

    snv_copynumber = parse_snv_cn.parse(somatic_calls, remixt)
    ideogram = read_ideogram.read()

    chromosomes = [chrom.lower() for chrom in chromosomes]

    for chrom in chromosomes:
        fig = plt.figure(constrained_layout=True, figsize=(15, 10))
        max = ideogram[ideogram.chrom == chrom].start.max()
        outer = fig.add_gridspec(1, 4, width_ratios=[max/250, 0.15, 0.15, 1-max/250])
        gs1 = gridspec.GridSpecFromSubplotSpec(9,1, height_ratios=[0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.25, 0.25], hspace=0.25,
                                               subplot_spec=outer[0])
        gs2 = gridspec.GridSpecFromSubplotSpec(9, 1, height_ratios=[0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.25, 0.25], hspace=0.25,
                                               subplot_spec=outer[1])
        gs3 = gridspec.GridSpecFromSubplotSpec(9, 1, height_ratios=[0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.25, 0.25], hspace=0.25,
                                               subplot_spec=outer[2])
        gs4 = gridspec.GridSpecFromSubplotSpec(9, 1, height_ratios=[0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.25, 0.25], hspace=0.25,
                                               subplot_spec=outer[3])

        axes = [plt.subplot(cell) for cell in gs1]
        axes.extend([plt.subplot(cell) for cell in gs2])
        axes.extend([plt.subplot(cell) for cell in gs3])
        axes.extend([plt.subplot(cell) for cell in gs4])
        
        for i, ax in enumerate(axes):
            if i not in [0,1,2,3,4,5,6,7,8,14]:
                ax.axis("off")

        axes = plot_chrom_on_axes(remixt, titan, roh, germline_calls,
                                  somatic_calls, tumour_coverage,
                                  normal_coverage, breakpoints, snv_copynumber,
                                  ideogram, chrom, axes)

        fig.suptitle("({}) Sample: {} Chromosome: {}".format("GRCh37", remixt_label, (chrom)))
        axes[8].set_xlabel("Position (Mb)", fontsize=14, fontname="Arial")

        axes = rasturize_axes(axes)

        plt.tight_layout()
        # write out
        pdf.savefig(fig)
    pdf.close()
