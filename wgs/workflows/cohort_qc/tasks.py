import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from wgs.utils import helpers
from classifycopynumber import parsers, transformations

def merge_segmental_cn(segmental_cn, concats):
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(concats, sep="\t", index=False)


def generate_segmental_copynumber(remixt, segmental_cn, sample):
    cn, stats = parsers.read_remixt(remixt)
    cn["sample"] = [sample] * len(cn)

    transformations.generate_segmental_cn(segmental_cn, cn, stats)


def merge_cna_tables(amps, dels, labels, output, cohort):
    cna = {d1: [a,d] for (d1, a), (d2, d) in zip(amps.items(), dels.items()) }
    number=0
    for s, files in cna.items():
        
        label = labels[(cohort, s)]
        for f in files:
        
            data = pd.read_csv(f, usecols=["cn_type", "gene_name", "pass_filter"])

            data=data[data.pass_filter == True]

            data=data.rename(columns={"cn_type": "CN", "gene_name":"Gene"})

            data["Sample_name"] = [label] * len(data)

            data = data[["Gene", "Sample_name", "CN"]]

            if number==0:
                header=True
            else:
                header=False

            if not data.empty:
                number+=1
                data.to_csv(output, index=False, mode='a', header=header, sep="\t")


def classify_remixt(sample_label, remixt, gtf, output_dir, amps, dels, docker_image=None):

    cmd = [
        "classifycopynumber", gtf, output_dir, sample_label, amps, dels, "--remixt_h5_filename", remixt, "--plot", False
    ]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)



def merge_mafs(mafs, merged_maf, write_header=True):
    for m in mafs:

        maf = pd.read_csv(m, sep="\t", chunksize=10e6)
        for chunk in maf:

            chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')
        #only write the first header
        write_header=False   


def merge_relabel_mafs(mafs, merged_maf, labels, class_label=None, write_header=True):

    mafs = list(mafs.values())
    labels = list(labels.values())

    #only write the first header
    #TODO better v names
    for label, maf in zip(labels, mafs):
        if maf:
            maf = pd.read_csv(maf, sep="\t", chunksize=10e6, skiprows=1)
            for chunk in maf:
                if labels:
                    chunk["Tumor_Sample_Barcode"] = [label] * len(chunk.index)
                if class_label:
                    chunk["Variant_Classification"] = chunk.Variant_Classification.apply(lambda c: c+class_label)
                chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')

            #only write the first header
            write_header=False
    

def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf, docker_image=None
):
    helpers.makedirs(tmpspace)

    cmd = [
        "/juno/work/shah/abramsd/oncokb-annotator/MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def label_variants_with_onco_annotatons(row):
    if row.oncogenic in ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]:
        return row.Variant_Classification +  "_" + row.oncogenic
    return row.Variant_Classification


def filter_annotated_maf(annotated_maf, filtered_maf, vcNames):
    annotated_maf = pd.read_csv(annotated_maf, sep="\t")
    annotated_maf["Variant_Classification"] = annotated_maf[["Variant_Classification", "oncogenic"]].apply(lambda row: 
        label_variants_with_onco_annotatons(row), axis=1
    )
    annotated_maf.to_csv(filtered_maf, sep="\t")

    classes = pd.DataFrame({"Variant_Classification":annotated_maf.Variant_Classification.unique()})
    classes.to_csv(vcNames, index=False)


def plot_mutation_burden(maf, burden_plot_path):
    maf = pd.read_csv(maf, sep="\t", usecols=["Tumor_Sample_Barcode"]).drop_duplicates()
    data = maf.groupby("Tumor_Sample_Barcode").size().sort_values(ascending=False)

    fig, axis = plt.subplots(figsize=(15, 5))
    nums = range(len(data.index))
    axis.bar(nums, data, align='center', color="black")
    plt.xticks(nums, data.index, rotation='vertical')
    axis.set_ylabel("Number of mutations")
    fig.savefig(burden_plot_path, format="png")
    plt.close()


def make_R_cohort_plots(
        cohort_maf, cntable, oncoplot_path, somatic_interactions,
        mafsummary, filtered_maf=None, vcNames=None, docker_image=None
):

    usemaf=filtered_maf
    if filtered_maf==None:
        usemaf=cohort_maf

    plots_cmd = [
        "maftools_plots.R", usemaf, vcNames, cntable,
        oncoplot_path, somatic_interactions, mafsummary
    ]

    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)


def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary, 
    burden_plot, report_path, docker_image=None
):
    cmd = [
        "run_cohort_qc_report.sh", report_path, cohort_label, oncoplot,
        somatic_interactions, mafsummary, burden_plot
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

