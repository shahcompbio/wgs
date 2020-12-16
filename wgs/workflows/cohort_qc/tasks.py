import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from wgs.utils import helpers
from classifycopynumber import parsers, transformations
import os


def merge_segmental_cn(segmental_cn, concats):
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(concats, sep="\t", index=False)


def generate_segmental_copynumber(remixt, segmental_cn, sample):

    cn, stats = parsers.read_remixt_parsed_csv(remixt)
    cn["sample"] = [sample] * len(cn)

    transformations.generate_segmental_cn(segmental_cn, cn, stats)


def merge_cna_tables(sample_labels, amps, dels, output, cohort):
    cna = {d1: [a,d] for (d1, a), (d2, d) in zip(amps.items(), dels.items()) }
    number=0
    for label, files in cna.items():
        matching_t_barcode = sample_labels[(cohort, label)]
        for f in files:
        
            data = pd.read_csv(f, usecols=["cn_type", "gene_name", "pass_filter"])

            data=data[data.pass_filter == True]

            data=data.rename(columns={"cn_type": "CN", "gene_name":"Gene"})

            data["Sample_name"] = matching_t_barcode

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
        "classifycopynumber", gtf, output_dir, sample_label, amps, dels, "--remixt_parsed_csv", remixt, "--plot", False
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def _write_maf(m, label, merged_maf, write_header):
    maf = pd.read_csv(m, sep="\t", chunksize=10e6)
    for chunk in maf:
        chunk["Tumor_Sample_Barcode"] = label
        chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a')
        write_header=False   


def merge_mafs(sample_labels, cohort, germline_mafs, somatic_mafs, merged_maf):
    write_header=True

    for label, m in germline_mafs.items():
        new_t_barcode = sample_labels[(cohort, label)]
        _write_maf(m, new_t_barcode, merged_maf, write_header)
        write_header=False

    for label, m in somatic_mafs.items():
        new_t_barcode = sample_labels[(cohort, label)]
        _write_maf(m, new_t_barcode, merged_maf, write_header)
        write_header=False


def annotate_maf(maf, annotated_maf, api_key, tempspace, class_label="", docker_image=None):

    annotate_maf_with_oncokb(maf, api_key, tmpspace, annotate_maf, docker_image)


def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf, docker_image=None
):
    helpers.makedirs(tmpspace)

    cmd = [
        "/juno/work/shah/abramsd/oncokb-annotator/MafAnnotator.py", "-i", maf, "-o", annotated_maf, "-b", api_key
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)



def filter_maf(annotated_maf, filtered_maf, write_header=True):
    oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
    maf = pd.read_csv(annotated_maf, sep="\t", chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.oncogenic.isin(oncogenic_annotations)]
        chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')

        write_header=False

def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
    maf = pd.read_csv(filtered_maf, sep="\t")
    maf["is_germline"] = [is_germline] * len(maf)
    maf.to_csv(annotated_maf, sep="\t", index=False)

def label_germline_somatic(row):
    if row.is_germline == True:
        return row.Variant_Classification + "_" + "germline"
    return row.Variant_Classification + "_" + "somatic"
  

def prepare_maf_for_maftools(cohort_label, filtered_maf, prepared_maf, non_synonymous_labels, vcNames):
    '''
    filter on non synonymous labels
    add germline/somatic annotate to variant classification
    add group/patient labels
    write out vcNames
    '''
    maf = pd.read_csv(filtered_maf, sep="\t")
    maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
    maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1 )
    nonsynclasses = pd.DataFrame({"Variant_Classification":maf.Variant_Classification.unique().tolist()})
    nonsynclasses.to_csv(vcNames, index=False)
    maf.to_csv(prepared_maf, sep="\t", index=False)


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
        maf, cntable, oncoplot_path, somatic_interactions,
        mafsummary, vcNames, docker_image=None
):

    plots_cmd = [
        "maftools_plots.R", maf, vcNames, cntable,
        oncoplot_path, somatic_interactions, mafsummary
    ]

    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)


def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary, 
    burden_plot, report_path, docker_image=None
):
    absolute_report = os.path.abspath(report_path)
    intermediate_dir = os.path.dirname(absolute_report)
    cmd = [
        "run_cohort_qc_report.sh", os.path.abspath(report_path), intermediate_dir, cohort_label, os.path.abspath(oncoplot),
        os.path.abspath(somatic_interactions), os.path.abspath(mafsummary), os.path.abspath(burden_plot)
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

