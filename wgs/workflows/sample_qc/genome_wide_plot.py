import matplotlib
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf

from wgs_qc_utils.plotter import variant_plotting
from wgs_qc_utils.plotter import coverage_plotting
from wgs_qc_utils.plotter import ideogram_plotting
from wgs_qc_utils.plotter import remixt_plotting
from wgs_qc_utils.plotter import roh_plotting
from wgs_qc_utils.plotter import snv_cn
from wgs_qc_utils.plotter import titan_plotting
from wgs_qc_utils.plotter import gene_annotation_plotting

from wgs_qc_utils.reader import read_remixt
from wgs_qc_utils.reader import read_roh
from wgs_qc_utils.reader import read_variant_calls
from wgs_qc_utils.reader import read_titan
from wgs_qc_utils.reader import read_coverage
from wgs_qc_utils.reader import parse_snv_cn

from wgs_qc_utils.reader.ideogram import read_ideogram




def plot_chrom_on_axes(
        remixt, titan, roh, germline_calls, somatic_calls,
        tumour_coverage, normal_coverage, breakpoints, vaf_data,
        ideogram, chrom, axes, normal_only=False
):
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
    prepped_roh = read_roh.prepare_at_chrom(roh, chrom)
    prepped_germline_calls = read_variant_calls.prepare_at_chrom(
        germline_calls, chrom, n_bins=2000
    )
    prepped_normal_coverage = read_coverage.prepare_at_chrom(normal_coverage, chrom)
    prepped_ideogram = read_ideogram.prepare_at_chrom(ideogram, chrom)

    coverage_ylim_min = 0
    coverage_ylim_max = 150
    coverage_cap_quantile = 0.9
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_max = 250
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_min = 0

    prepped_normal_coverage_cap = prepped_normal_coverage.coverage.quantile(
            coverage_cap_quantile
        )

    if not normal_only:
        prepped_remxit = read_remixt.prepare_at_chrom(remixt, chrom)
        prepped_snv_cn = parse_snv_cn.prepare_at_chrom(vaf_data, chrom)
        prepped_titan = read_titan.prepare_at_chrom(titan, chrom)
        prepped_somatic_calls = read_variant_calls.prepare_at_chrom(
            somatic_calls, chrom, n_bins=2000
        )
        prepped_tumour_coverage = read_coverage.prepare_at_chrom(tumour_coverage, chrom)
        prepped_breakpoints = read_variant_calls.prepare_at_chrom(
            breakpoints, chrom, n_bins=2000
        )
        prepped_tumour_coverage_cap = prepped_tumour_coverage.coverage.quantile(
                coverage_cap_quantile
        )
        if prepped_tumour_coverage_cap > coverage_ylim_max:
            coverage_ylim_max = 50 * ((int(prepped_tumour_coverage_cap) / 50) + 1)

    if prepped_normal_coverage_cap > coverage_ylim_max:
        coverage_ylim_max = 50 * ((int(prepped_normal_coverage_cap) / 50) + 1)

    anno_genes = gene_annotation_plotting.get_gene_annotation_data(chrom)
    chrom_max = prepped_ideogram.start.max()
    if normal_only:

        axes[0] = variant_plotting.plot_fill(
            prepped_germline_calls.location, prepped_germline_calls.n_events,
            axes[0], "Germline \n Alteration \n Frequency", chrom_max, "germline"
        )
        axes[1] = coverage_plotting.plot(
            prepped_normal_coverage.start, prepped_normal_coverage.coverage,
            coverage_ylim_min, coverage_ylim_max,
            axes[1], "Normal \n Coverage", chrom_max, anno_genes=anno_genes
        )
        axes[2] = roh_plotting.plot(
            prepped_roh.start, prepped_roh.state, axes[2], chrom_max)

        axes[3] = ideogram_plotting.plot(prepped_ideogram, axes[3])

        if not anno_genes.empty:
            axes[5] = gene_annotation_plotting.add_gene_annotation_legend(anno_genes, axes[5])

    else:
        axes[0] = variant_plotting.plot_fill(
            prepped_somatic_calls.location, prepped_somatic_calls.n_events,
            axes[0], "Somatic \n Alteration \n Frequency", chrom_max, "somatic"
        )
        axes[1] = variant_plotting.plot_fill(
            prepped_germline_calls.location, prepped_germline_calls.n_events,
            axes[1], "Germline \n Alteration \n Frequency", chrom_max, "germline"
        )
        axes[2] = coverage_plotting.plot(
            prepped_tumour_coverage.start, prepped_tumour_coverage.coverage,
            coverage_ylim_min, coverage_ylim_max,
            axes[2], "Tumour \n Coverage", chrom_max
        )
        axes[3] = coverage_plotting.plot(
            prepped_normal_coverage.start, prepped_normal_coverage.coverage,
            coverage_ylim_min, coverage_ylim_max,
            axes[3], "Normal \n Coverage", chrom_max
        )

        axes[4] = variant_plotting.plot_bar(
            prepped_breakpoints.location, prepped_breakpoints.n_events,
            axes[4], "Breakpoint \n Frequency", chrom_max
        )

        axes[5] = remixt_plotting.plot(
            prepped_remxit.start, prepped_remxit.major_raw,
            prepped_remxit.minor_raw, axes[5], chrom_max
        )
        axes[5] = snv_cn.plot_scatter(
            prepped_snv_cn.pos, prepped_snv_cn.frac_cn, axes[5])

        axes[6] = titan_plotting.plot(
            prepped_titan.Position, prepped_titan.LogRatio,
            prepped_titan.color, axes[6], chrom_max, anno_genes=anno_genes)

        axes[7] = roh_plotting.plot(
            prepped_roh.start, prepped_roh.state, axes[7], chrom_max)

        axes[8] = ideogram_plotting.plot(prepped_ideogram, axes[8])

        axes[14] = snv_cn.plot_hist(prepped_snv_cn.frac_cn, axes[14])

        axes[23] = remixt_plotting.add_remixt_legend(axes[23])
        axes[15] = titan_plotting.add_titan_legend(axes[15])

        if not anno_genes.empty:
            axes[24] = gene_annotation_plotting.add_gene_annotation_legend(anno_genes, axes[24])

    return axes


def rasterize_axes(axes):
    """
    rasterize axes
    :param axes: list of axes to rasterize
    :return: rasterized axes
    """
    for ax in axes:
        ax.set_rasterized(True)
    return axes


def _add_axis_unit(frame, normal_only=False):
    if normal_only:
        return gridspec.GridSpecFromSubplotSpec(
            4, 1, height_ratios=[0.5, 0.5, 0.25, 0.1],
            hspace=0.25,
            subplot_spec=frame
        )
    return gridspec.GridSpecFromSubplotSpec(
        9, 1, height_ratios=[0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.25, 0.25],
        hspace=0.25,
        subplot_spec=frame
    )

def _make_axes(ideogram, chrom, sample, fig, normal_only=False):
    max = ideogram[ideogram.chrom == chrom].start.max()
    outer = fig.add_gridspec(1, 4, width_ratios=[max / 250, 0.15, 0.15, 1 - max / 250])

    gs1 = _add_axis_unit(outer[0], normal_only=normal_only)
    gs2 = _add_axis_unit(outer[1], normal_only=normal_only)
    gs3 = _add_axis_unit(outer[2], normal_only=normal_only)
    gs4 = _add_axis_unit(outer[3], normal_only=normal_only)

    axes = [plt.subplot(cell) for cell in gs1]
    axes.extend([plt.subplot(cell) for cell in gs2])
    axes.extend([plt.subplot(cell) for cell in gs3])
    axes.extend([plt.subplot(cell) for cell in gs4])

    for i, ax in enumerate(axes):
        if normal_only:
            if i not in [0, 1, 2, 3]:
                ax.axis("off")
        else:
            if i not in [0, 1, 2, 3, 4, 5, 6, 7, 8, 14]:
                ax.axis("off")

    fig.suptitle("({}) Sample: {} Chromosome: {}".format("GRCh37", sample, (chrom)))
    if normal_only:
        axes[3].set_xlabel("Position (Mb)", fontsize=14, fontname="Arial")
    else:
        axes[8].set_xlabel("Position (Mb)", fontsize=14, fontname="Arial")
    return axes


def genome_wide_plot(
        remixt, remixt_label, titan, roh, germline_calls, somatic_calls,
        tumour_coverage, normal_coverage, breakpoints, chromosomes, pdf,
        normal_only=False
):
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
    roh = read_roh.read(roh)
    germline_calls = read_variant_calls.read(germline_calls)
    normal_coverage = read_coverage.read(normal_coverage)
    if normal_only:
        remixt = None
        somatic_calls = None
        titan = None
        tumour_coverage = None
        breakpoints =None
        snv_copynumber = None
    else:
        remixt = read_remixt.read(remixt)
        somatic_calls = read_variant_calls.read(somatic_calls)
        titan = read_titan.read(titan)
        tumour_coverage = read_coverage.read(tumour_coverage)
        breakpoints = read_variant_calls.read_svs(breakpoints)
        snv_copynumber = parse_snv_cn.parse(somatic_calls, remixt)

    ideogram = read_ideogram.read()

    chromosomes = [chrom.lower() for chrom in chromosomes]

    for chrom in chromosomes:
        fig = plt.figure(constrained_layout=True, figsize=(15, 10))

        axes = _make_axes(ideogram, chrom, remixt_label, fig, normal_only=normal_only)

        axes = plot_chrom_on_axes(remixt, titan, roh, germline_calls,
                                  somatic_calls, tumour_coverage,
                                  normal_coverage, breakpoints, snv_copynumber,
                                  ideogram, chrom, axes, normal_only=normal_only)

        rasterize_axes(axes)

        plt.tight_layout()
        pdf.savefig(fig)
    pdf.close()
