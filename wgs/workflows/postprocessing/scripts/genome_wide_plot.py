import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf

from . import roh_plotting
from . import variant_plotting
from . import coverage_plotting
from . import titan_plotting
from . import ideogram_plotting
from . import snv_cn
from . import remixt_plotting
from . import gene_annotation_plotting


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
    prepped_remxit = remixt_plotting.prepare_at_chrom(remixt, chrom)
    prepped_titan = titan_plotting.prepare_at_chrom(titan, chrom)
    prepped_roh = roh_plotting.prepare_at_chrom(roh, chrom)
    prepped_somatic_calls = variant_plotting.prepare_at_chrom(somatic_calls, chrom, n_bins=2000)
    prepped_germline_calls = variant_plotting.prepare_at_chrom(germline_calls, chrom, n_bins=2000)
    prepped_tumour_coverage = coverage_plotting.prepare_at_chrom(tumour_coverage, chrom, bin=True, n_bins=2000)
    prepped_normal_coverage = coverage_plotting.prepare_at_chrom(normal_coverage, chrom, bin=True, n_bins=2000)
    prepped_breakpoints = variant_plotting.prepare_at_chrom(breakpoints, chrom, n_bins=2000)
    prepped_vaf = snv_cn.prepare_at_chrom(vaf_data, chrom)
    prepped_ideogram = ideogram_plotting.prepare_at_chrom(ideogram, chrom)

    coverage_ylim_max = prepped_tumour_coverage.coverage.max() + 10
    coverage_ylim_min = prepped_normal_coverage.coverage.min() - 10
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_max = 250
    if pd.isnull(coverage_ylim_max):
        coverage_ylim_min = 250
    anno_genes = gene_annotation_plotting.get_gene_annotation_data(chrom)

    chrom_max = prepped_ideogram.start.max()
    axes[0] = variant_plotting.plot_fill(prepped_somatic_calls, axes[0], "Somatic \n Alteration \n Frequency", chrom_max)
    axes[1] = variant_plotting.plot_fill(prepped_germline_calls, axes[1], "Germline \n Alteration \n Frequency", chrom_max)
    axes[2] = coverage_plotting.plot(prepped_tumour_coverage,
                                     coverage_ylim_min, coverage_ylim_max,
                                     axes[2], "Tumour \n Coverage", chrom_max)
    axes[3] = coverage_plotting.plot(prepped_normal_coverage,
                                     coverage_ylim_min, coverage_ylim_max,
                                     axes[3], "Normal \n Coverage", chrom_max)
    axes[4] = variant_plotting.plot_bar(prepped_breakpoints, axes[4], "Breakpoint \n Frequency", chrom_max)


    axes[5] = remixt_plotting.plot(prepped_remxit, axes[5], chrom_max)
    axes[5] = snv_cn.plot_scatter(prepped_vaf, axes[5])

    axes[6] = titan_plotting.plot(prepped_titan, anno_genes, axes[6], chrom_max)
    axes[7] = roh_plotting.plot(prepped_roh, axes[7])

    axes[8] = ideogram_plotting.plot(prepped_ideogram, axes[8])

    axes[14] = snv_cn.plot_hist(prepped_vaf, axes[14])

    axes[23] = remixt_plotting.add_remixt_legend(axes[23])
    axes[15] = titan_plotting.add_titan_legend( axes[15])

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
                     tumour_coverage, normal_coverage, breakpoints, ideogram, chromosomes, pdf):

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
    remixt = remixt_plotting.read(remixt, remixt_label)
    titan = titan_plotting.read(titan)
    roh = roh_plotting.read(roh)
    germline_calls = variant_plotting.read(germline_calls)
    tumour_coverage = coverage_plotting.read(tumour_coverage)

    normal_coverage = coverage_plotting.read(normal_coverage)
    breakpoints = pd.read_csv(breakpoints)[["chromosome_1", "chromosome_2", "position_1", "position_2"]]
    breakpoints = breakpoints.astype({"chromosome_1": str, "chromosome_2": str})

    breakpoints = pd.DataFrame({"chr": breakpoints["chromosome_1"].append(breakpoints["chromosome_2"]),
                                "pos": breakpoints["position_1"].append(breakpoints["position_2"])})

    somatic_calls = variant_plotting.read_consensus(somatic_calls)

    vaf_data = snv_cn.parse(somatic_calls, remixt)
    ideogram = ideogram_plotting.read(ideogram)
    # get chroms to plot on and make sure all data contains those chroms

    # for chrom in map(str, list(range(1,23)) + ["X"]):
    for chrom in chromosomes:
        fig = plt.figure(constrained_layout=True, figsize=(15, 10))

        ref = 250
        MAX = ideogram[ideogram.chrom == chrom].start.max()
        perc_off = MAX
        outer = fig.add_gridspec(1, 4, width_ratios=[MAX/ref, 0.15, 0.15, 1-MAX/ref])
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
        #
        for i, ax in enumerate(axes):
            if i not in [0,1,2,3,4,5,6,7,8,14]:
                ax.axis("off")

        axes = plot_chrom_on_axes(remixt, titan, roh, germline_calls,
                                  somatic_calls, tumour_coverage,
                                  normal_coverage, breakpoints, vaf_data,
                                  ideogram, chrom, axes)

        fig.suptitle("({}) Sample: {} Chromosome: {}".format("GRCh37", remixt_label, (chrom)))
        axes[8].set_xlabel("Position (Mb)", fontsize=14, fontname="Arial")

        axes = rasturize_axes(axes)

        plt.tight_layout()
        # write out
        print(pdf)
        pdf.savefig(fig)
    pdf.close()
