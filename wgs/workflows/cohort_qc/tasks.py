import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from wgs.utils import helpers


def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf, docker_image=None
):
    helpers.makedirs(tmpspace)

    cmd = [
        "MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def filter_annotated_maf(annotated_maf, filtered_maf):
    annotated_maf = pd.read_csv(annotated_maf, sep="\t")

    filt_maf = annotated_maf[
        (annotated_maf.ONCOGENIC == "Oncogenic")
        | (annotated_maf.ONCOGENIC == "Likely Oncogenic")
        ]
    filt_maf.to_csv(filtered_maf, sep="\t")


def plot_mutation_burden(maf, burden_plot_path):
    maf = pd.read_csv(maf, sep="\t").drop_duplicates()
    data = maf.groupby("Tumor_Sample_Barcode").size().sort_values(ascending=False)

    fig, axis = plt.subplots(figsize=(15, 5))
    nums = range(len(data.index))
    axis.bar(nums, data, align='center', color="black")
    plt.xticks(nums, data.index, rotation='vertical')
    axis.set_ylabel("Number of mutations")
    fig.savefig(burden_plot_path, format="png")
    plt.close()


def make_R_cohort_plots(
        cohort_maf, oncoplot_path, somatic_interactions,
        mafsummary, filtered_maf=None
):
    if not filtered_maf:
        filtered_maf = cohort_maf

    plots_cmd = [
        "maftools_plots.R", cohort_maf, filtered_maf,
        oncoplot_path, somatic_interactions, mafsummary
    ]

    pypeliner.commandline.execute(*plots_cmd)


def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary, burden_plot, report_path):
    cmd = [
        "run_report.sh", report_path, cohort_label, oncoplot,
        somatic_interactions, mafsummary, burden_plot
    ]
    pypeliner.commandline.execute(*cmd)
