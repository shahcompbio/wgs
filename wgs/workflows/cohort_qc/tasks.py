import pandas as pd
import pypeliner
from wgs.utils import helpers
import matplotlib.pyplot as plt

# def group_maf(maf, groups, groupped_maf, label_col="Tumor_Sample_Barcode"):
#     maf = pd.read_csv(maf, sep="\t")
#     groups = pd.read_csv(groups).set_index("label")
#     maf[label_col] = maf.replace(groups.to_dict())
#     maf.to_csv(groupped_maf, sep="\t", index=False)
    
def annotate_maf_with_oncokb(maf, api_key, tmpspace, annotated_maf, docker_image=None):
    helpers.makedirs(tmpspace)

    cmd = ["python", "/juno/work/shah/abramsd/oncokb-annotator/MafAnnotator.py", "-i",
        maf, "-o", annotated_maf, "-b", api_key]#"86c693fc-7cac-4032-9c32-45e07f7e0965 "] ####!!!!TODO: remove API key

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def filter_annotated_maf(annotated_maf, filtered_maf):
    
    annotated_maf = pd.read_csv(annotated_maf, sep="\t")

    filt_maf = annotated_maf[(annotated_maf.oncogenic == "Oncogenic") 
        | (annotated_maf.oncogenic == "Likely Oncogenic")]
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


def make_R_cohort_plots(cohort_maf, oncoplot_path, somatic_interactions, 
    mafsummary, filtered_maf=None
):

    if not filtered_maf:
        filtered_maf = cohort_maf

    #TODO: add to conda env - had issues wiht adding maftools
    plots_cmd = ["Rscript", "/juno/work/shah/abramsd/CODE/wgs/wgs/workflows/cohort_qc/scripts/plots.R", 
        cohort_maf, filtered_maf, oncoplot_path, somatic_interactions, mafsummary]

    pypeliner.commandline.execute(*plots_cmd)


def make_report(cohort_label, oncoplot,somatic_interactions,mafsummary, burden_plot, report_path):
    
    # cmd = [ "R", "-e", "\"", "library(rmarkdown)", ";", "render('/juno/work/shah/abramsd/CODE/wgs/wgs/workflows/cohort_qc/scripts/report.Rmd', output_file='{}', params = list(label='{}', oncoplot='{}', somatic_plot='{}', summary='{}'))\"".format(report_path, cohort_label, oncoplot, somatic_interactions, mafsummary)]
    cmd = ["run_report.sh", report_path, cohort_label, oncoplot, 
    somatic_interactions, mafsummary, burden_plot]
    pypeliner.commandline.execute(*cmd)
