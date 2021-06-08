# import matplotlib.pyplot as plt
# import pandas as pd
# import pypeliner
# from wgs.utils import helpers
# from classifycopynumber import parsers, transformations
# import os
# import wgs_analysis.algorithms.cnv
# import mafannotator.MafAnnotator as ma
# import numpy as np
# from rpy2.robjects.packages import importr
# import rpy2.robjects as robjects
# rmarkdown = importr("rmarkdown")
# import seaborn
#
#
# def build_gene_list(cna_table, genelist, base_genes='default'):
#     if base_genes=='default':
#         genes=["PPM1D", "TP53",  "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2",
#         "ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC"]
#     genes = pd.DataFrame({"gene":list(genes)})
#     genes.to_csv(genelist)
#
#
# def merge_segmental_cn(segmental_cn, concats):
#     files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
#     segmental_cn_combined = pd.concat(files)
#     segmental_cn_combined.to_csv(concats, sep="\t", index=False)
#


#
# def generate_segmental_copynumber(remixt, segmental_cn, sample):
#
#     cn, stats = parsers.read_remixt_parsed_csv(remixt)
#     cn = cn.astype({"chromosome":"str"})
#
#     cn["sample"] = [sample] * len(cn)
#
#     aggregated_cn_data = wgs_analysis.algorithms.cnv.aggregate_adjacent(
#         cn,
#         value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
#         stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2', 'sample'],
#         length_normalized_cols=['major_raw', 'minor_raw'],
#     )
#     aggregated_cn_data['copy'] = aggregated_cn_data['major_raw'] + aggregated_cn_data['minor_raw']
#
#
#     generate_segmental_cn(segmental_cn, aggregated_cn_data, stats["ploidy"])
#
#
# def generate_gistic_outputs(gistic_data, hdel_data, cbio_table):
#     gistic_data['gistic_value'] = 2
#     gistic_data.loc[gistic_data['log_change'] < 1, 'gistic_value'] = 1
#     gistic_data.loc[gistic_data['log_change'] < 0.5, 'gistic_value'] = 0
#     gistic_data.loc[gistic_data['log_change'] < -0.5, 'gistic_value'] = -1
#
#     # Merge hdels
#     hdel_data['is_hdel'] = 1
#     gistic_data = gistic_data.merge(hdel_data[['Hugo_Symbol', 'sample', 'is_hdel']], how='left')
#     gistic_data['is_hdel'] = gistic_data['is_hdel'].fillna(0).astype(int)
#     gistic_data.loc[gistic_data['is_hdel'] == 1, 'gistic_value'] = -2
#
#     # Gistic_data generation
#     gistic_data = gistic_data[['Hugo_Symbol', 'sample', 'gistic_value']]
#     gistic_data = gistic_data.drop_duplicates()
#     gistic_data = gistic_data.astype({"gistic_value": "Int64"})
#
#     gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
#     gistic_matrix.reset_index(inplace=True)
#     gistic_matrix.to_csv(cbio_table, sep="\t", index=False, na_rep="NA")
#
#
# def make_cbio_cna_table(amps, dels, cbio_table):
#     amps = pd.read_csv(amps,  sep="\t", usecols=["gene_name", "log_change", "sample"])
#     amps = amps.rename(columns={"gene_name":"Hugo_Symbol"})
#
#     dels = pd.read_csv(dels,  sep="\t", usecols=["gene_name", "sample"])
#     dels = dels.rename(columns={"gene_name":"Hugo_Symbol"})
#
#     generate_gistic_outputs(amps, dels, cbio_table)
#
#
# def make_maftools_cna_table(amps, dels, maftools_table):
#     print(amps)
#     amps = pd.read_csv(amps, sep="\t", usecols=["gene_name", "sample", "pass_filter"])
#     amps = amps.rename(columns={"gene_name":"Gene", "sample": "Sample_name"})
#     amps["CN"] = "amplification"
#     amps=amps[amps.pass_filter == True]
#
#     dels = pd.read_csv(dels,  sep="\t", usecols=["gene_name", "sample", "pass_filter"])
#     dels = dels.rename(columns={"gene_name":"Gene", "sample": "Sample_name"})
#     dels["CN"] = "deletion"
#     dels=dels[dels.pass_filter == True]
#
#     out = pd.concat([amps, dels])
#     out = out[["Gene",  "Sample_name", "CN"]]
#     out.to_csv(maftools_table, index=False, sep="\t")
#
#
# def merge_cna_tables(tables, output):
#     number=0
#     for label, cna in tables.items():
#
#         data = pd.read_csv(cna)
#         data["sample"] = label
#         if number==0:
#             header=True
#         else:
#             header=False
#         number+=1
#         data.to_csv(output, index=False, mode='a', header=header, sep="\t")
#
#
#
#
# def classify_remixt(sample_label, remixt, gtf, output_dir, amps, dels, docker_image=None):
#
#     cmd = [
#         "classifycopynumber", gtf, output_dir, amps, dels, "--remixt_parsed_csv", remixt
#     ]
#     pypeliner.commandline.execute(*cmd, docker_image=docker_image)
#
#
# def _write_maf(m, label, merged_maf, write_header):
#     maf = pd.read_csv(m, sep="\t", dtype='str', chunksize=10e6)
#     print(m)
#     for chunk in maf:
#         chunk["Tumor_Sample_Barcode"] = label
#         print(chunk)
#         print(chunk[["t_ref_count", "t_alt_count",
#                                  "n_ref_count", "n_alt_count"]].dtypes)
#         chunk["t_ref_count"] = pd.to_numeric(chunk['t_ref_count'],
#             downcast='float', errors='raise').astype('Int64')
#         chunk["t_alt_count"] = pd.to_numeric(chunk['t_alt_count'],
#             downcast='float', errors='raise').astype('Int64')
#         chunk["n_ref_count"] = pd.to_numeric(chunk['n_ref_count'],
#             downcast='float', errors='raise').astype('Int64')
#         chunk["n_alt_count"] = pd.to_numeric(chunk['n_alt_count'],
#             downcast='float', errors='raise').astype('Int64')
#
#         print(chunk[["t_ref_count", "t_alt_count",
#                                  "n_ref_count", "n_alt_count"]].dtypes)
#
#         #chunk = chunk.astype({"t_ref_count":"Int64", "t_alt_count":"Int64",
#          #   "n_ref_count":"Int64", "n_alt_count":"Int64"})
#
#         chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a', na_rep="")
#         write_header=False
#
#
# def merge_mafs(germline, somatic_mafs, merged_maf):
#     write_header=True
#
#     for label, m in somatic_mafs.items():
#         _write_maf(m, label, merged_maf, write_header)
#         write_header=False
#     for label, m in germline.items():
#         _write_maf(m, label, merged_maf, write_header)
#         write_header=False
#
# def annotate_maf(maf, annotated_maf, api_key, tempspace, class_label="", docker_image=None):
#
#     annotate_maf_with_oncokb(maf, api_key, tmpspace, annotate_maf, docker_image)
#
#
# def annotate_maf_with_oncokb(
#         maf, api_key, tmpspace, annotated_maf, docker_image=None
# ):
#     helpers.makedirs(tmpspace)
#     ma.annotate(maf, annotated_maf, api_key)
#
#
# def filter_maf(annotated_maf, filtered_maf, write_header=True):
#     oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
#     maf = pd.read_csv(annotated_maf, sep="\t", dtype='str', chunksize=10e6)
#     for chunk in maf:
#
#         chunk = chunk[chunk.ONCOGENIC.isin(oncogenic_annotations)]
#         chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')
#
#         write_header=False
#
#
# def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
#     maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
#     maf["is_germline"] = [is_germline] * len(maf)
#     maf.to_csv(annotated_maf, sep="\t", index=False)
#
#
# def label_germline_somatic(row):
#     if row.is_germline == True:
#         return row.Variant_Classification + "_" + "germline"
#     return row.Variant_Classification + "_" + "somatic"
#
#
# def prepare_maf_for_maftools(cohort_label, filtered_maf, prepared_maf, non_synonymous_labels, vcNames):
#     '''
#     filter on non synonymous labels
#     add germline/somatic annotate to variant classification
#     add group/patient labels
#     write out vcNames
#     '''
#     maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
#     maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
#     maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1 )
#     nonsynclasses = pd.DataFrame({"Variant_Classification":maf.Variant_Classification.unique().tolist()})
#     nonsynclasses.to_csv(vcNames, index=False)
#     maf.to_csv(prepared_maf, sep="\t", index=False)
#
#
# def plot_mutation_burden(maf, burden_plot_path):
#     maf = pd.read_csv(maf, dtype='str', sep="\t").drop_duplicates()
#     key_genes= ["PPM1D", "TP53",  "BRCA1", "BRCA2", "MECOM",
#         "RB1", "PTEN", "PALB2","ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC"]
#
#     key_genes_maf = maf[maf.Hugo_Symbol.isin(key_genes)]
#
#     main_counts = maf.groupby("Tumor_Sample_Barcode").size()
#     key_genes_counts = key_genes_maf.groupby("Tumor_Sample_Barcode").size()
#
#     fig, axis = plt.subplots(1, 2, figsize=(20, 10))
#
#     seaborn.barplot(x=main_counts.index.tolist(), y=main_counts.tolist(), ax=axis[0])
#     axis[0].set_ylabel("Number of mutations")
#     plt.sca(axis[0])
#     plt.xticks(rotation=90)
#     axis[0].tick_params(axis='both', which='major', labelsize=8)
#     axis[0].set_title("All genes")
#
#     seaborn.barplot(x=key_genes_counts.index.tolist(), y=key_genes_counts.tolist(), ax=axis[1])
#     axis[1].set_ylabel("Number of mutations")
#     plt.sca(axis[1])
#     plt.xticks(rotation=90)
#     axis[1].tick_params(axis='both', which='major', labelsize=8)
#     axis[1].set_title("key genes: {}".format(" ".join(key_genes)))
#
#     fig.subplots_adjust( bottom=0.45)
#     fig.savefig(burden_plot_path, format="png")
#
#     plt.close()
#
# def make_R_cohort_plots(
#         maf, cntable, oncoplot_path, somatic_interactions,
#         mafsummary, vcNames, genelist, docker_image=None
# ):
#     maftools_plots_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
#         'scripts','maftools_plots.R'
#     )
#     plots_cmd = [
#         "Rscript", maftools_plots_script, maf, vcNames, cntable,
#         oncoplot_path, somatic_interactions, mafsummary, genelist
#     ]
#
#     pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)
#
#
# def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary,
#     burden_plot, report_path, docker_image=None
# ):
#     rmd_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
#         'scripts','report.Rmd'
#     )
#
#     parameters = robjects.r.list(label=cohort_label,
#         oncoplot=os.path.abspath(oncoplot), somatic_plot= os.path.abspath(somatic_interactions),
#         summary=os.path.abspath(mafsummary), burden_plot=os.path.abspath(burden_plot)
#     )
#
#     rmarkdown.render(rmd_script, output_file= os.path.abspath(report_path),
#         intermediates_dir = os.path.dirname(os.path.abspath(report_path)),
#         output_options=robjects.r.list(self_contained=True), params=parameters
#     )
#
#
#     #cmd = [
#     #    "run_cohort_qc_report.sh", os.path.abspath(report_path), intermediate_dir, cohort_label, os.path.abspath(oncoplot),
#     #    os.path.abspath(somatic_interactions), os.path.abspath(mafsummary), os.path.abspath(burden_plot)
#     #]
#     #pypeliner.commandline.execute(*cmd, docker_image=docker_image)
#

import matplotlib.pyplot as plt
import pandas as pd
import pypeliner
from wgs.utils import helpers
from classifycopynumber import parsers, transformations
import os
import wgs_analysis.algorithms.cnv
import mafannotator.MafAnnotator as ma
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
rmarkdown = importr("rmarkdown")

def build_gene_list(cna_table, genelist, base_genes='default'):
    if base_genes == 'default':
        genes = ["PPM1D", "TP53", "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2",
                 "ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC"]
    genes = pd.DataFrame({"gene": list(genes)})
    genes.to_csv(genelist)


def merge_segmental_cn(segmental_cn, concats):
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(concats, sep="\t", index=False)

def generate_segmental_cn(filename, aggregated_cn_data, ploidy,  cn_col="copy", length_col="length"):
    aggregated_cn_data['ploidy'] = ploidy
    aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data[cn_col] / aggregated_cn_data['ploidy'])
    aggregated_cn_data['num.mark'] = (aggregated_cn_data[length_col] / 500000).astype(int)
    aggregated_cn_data = aggregated_cn_data.rename(columns={'sample': 'ID', 'chromosome': 'chrom', 'start': 'loc.start', 'end': 'loc.end'})
    aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
    aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
    aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
    aggregated_cn_data = transformations._correct_seg_bin_ends(aggregated_cn_data)
    aggregated_cn_data.to_csv(filename, index=None, sep='\t')

def generate_segmental_copynumber(remixt, segmental_cn, sample):
    cn, stats = parsers.read_remixt_parsed_csv(remixt)
    cn = cn.astype({"chromosome": "str"})

    cn["sample"] = [sample] * len(cn)

    aggregated_cn_data = wgs_analysis.algorithms.cnv.aggregate_adjacent(
        cn,
        value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
        stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2', 'sample'],
        length_normalized_cols=['major_raw', 'minor_raw'],
    )
    aggregated_cn_data['copy'] = aggregated_cn_data['major_raw'] + aggregated_cn_data['minor_raw']

    generate_segmental_cn(segmental_cn, aggregated_cn_data, stats["ploidy"])


def make_cbio_cna_table(cn_change_filename, cbio_table):
    gistic_data = pd.read_csv(cn_change_filename, sep="\t", usecols=["gene_name", "sample", "gistic_value"])
    gistic_data = gistic_data.rename(columns={"gene_name": "Hugo_Symbol"})
    gistic_data = gistic_data.drop_duplicates()
    gistic_data = gistic_data.astype({"gistic_value": "Int64"})

    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)
    gistic_matrix.to_csv(cbio_table, sep="\t", index=False, na_rep="NA")


def make_maftools_cna_table(cn_change_filename, maftools_table):
    cn_change = pd.read_csv(cn_change_filename, sep="\t",
                            usecols=["gene_name", "sample", "is_hdel", "is_loh", "is_hlamp"])

    cn_change = cn_change[cn_change['is_hdel'] | cn_change['is_loh'] | cn_change['is_hlamp']]

    cn_change["is_loh"] = cn_change.is_loh.replace({False: "", True: "loh"})
    cn_change["is_hdel"] = cn_change.is_hdel.replace({False: "", True: "hdel"})
    cn_change["is_hlamp"] = cn_change.is_hlamp.replace({False: "", True: "hlamp"})
    cn_change["CN"] = cn_change[["is_hdel", "is_loh", "is_hlamp"]].apply(
        lambda row: "".join(map(str, row)) ,axis=1
    )

    cn_change = cn_change.rename(columns={"gene_name": "Gene", "cn_type": "CN", "sample": "Sample_Name"})
    cn_change = cn_change[["Gene", "Sample_Name", "CN"]]
    cn_change.to_csv(maftools_table, index=False, sep="\t")


def merge_cna_tables(tables, output):
    number = 0
    for label, cna in tables.items():

        data = pd.read_csv(cna)
        data["sample"] = label
        if number == 0:
            header = True
        else:
            header = False
        number += 1
        data.to_csv(output, index=False, mode='a', header=header, sep="\t")


def classify_remixt(sample_label, remixt, gtf, cn_change, docker_image=None):
    cmd = [
        "classifycopynumber", gtf, cn_change, "--remixt_parsed_csv", remixt, '--sample_ids', sample_label,
    ]
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


# def _write_maf(m, label, merged_maf, write_header):
#     maf = pd.read_csv(m, sep="\t", dtype='str', chunksize=10e6)
#     for chunk in maf:
#         chunk["Tumor_Sample_Barcode"] = label
#         print(chunk["t_ref_count"], chunk.t_ref_count.unique())
#         print(chunk.t_ref_count.unique())
#         chunk.t_ref_count.astype({"t_ref_count": "int64"})
#         cc
#         # chunk = chunk.astype({"t_ref_count": pd.Int64Dtype()})
#         # , "t_alt_count": pd.Int64Dtype(),
#         #                       "n_ref_count": pd.Int64Dtype(), "n_alt_count": pd.Int64Dtype()})
#
#         chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a', na_rep="")
#         write_header = False
def _write_maf(m, label, merged_maf, write_header):
    maf = pd.read_csv(m, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk["Tumor_Sample_Barcode"] = label
        chunk["t_ref_count"] = pd.to_numeric(chunk['t_ref_count'],
            downcast='float', errors='raise').astype('Int64')
        chunk["t_alt_count"] = pd.to_numeric(chunk['t_alt_count'],
            downcast='float', errors='raise').astype('Int64')
        chunk["n_ref_count"] = pd.to_numeric(chunk['n_ref_count'],
            downcast='float', errors='raise').astype('Int64')
        chunk["n_alt_count"] = pd.to_numeric(chunk['n_alt_count'],
            downcast='float', errors='raise').astype('Int64')

        chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, mode='a', na_rep="")
        write_header=False


def merge_mafs(germline, somatic_mafs, merged_maf):
    write_header = True

    for label, m in somatic_mafs.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header = False
    for label, m in germline.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header = False


def annotate_maf(maf, annotated_maf, api_key, tempspace, class_label="", docker_image=None):
    annotate_maf_with_oncokb(maf, api_key, tmpspace, annotate_maf, docker_image)


def annotate_maf_with_oncokb(
        maf, api_key, tmpspace, annotated_maf, docker_image=None
):
    helpers.makedirs(tmpspace)
    ma.annotate(maf, annotated_maf, api_key)


def filter_maf(annotated_maf, filtered_maf, write_header=True):
    oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
    maf = pd.read_csv(annotated_maf, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.ONCOGENIC.isin(oncogenic_annotations)]
        chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')

        write_header = False


def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
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
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
    maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
    maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1)
    nonsynclasses = pd.DataFrame({"Variant_Classification": maf.Variant_Classification.unique().tolist()})
    nonsynclasses.to_csv(vcNames, index=False)
    maf.to_csv(prepared_maf, sep="\t", index=False)


def plot_mutation_burden(maf, burden_plot_path):
    maf = pd.read_csv(maf, dtype='str', sep="\t", usecols=["Tumor_Sample_Barcode"]).drop_duplicates()
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
        mafsummary, vcNames, genelist, docker_image=None
):
    maftools_plots_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'scripts','maftools_plots.R'
    )
    plots_cmd = [
        "Rscript", maftools_plots_script, maf, vcNames, cntable,
        oncoplot_path, somatic_interactions, mafsummary, genelist
    ]

    pypeliner.commandline.execute(*plots_cmd, docker_image=docker_image)



def make_report(cohort_label, oncoplot, somatic_interactions, mafsummary,
    burden_plot, report_path, docker_image=None
):
    rmd_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'scripts','report.Rmd'
    )

    parameters = robjects.r.list(label=cohort_label,
        oncoplot=os.path.abspath(oncoplot), somatic_plot= os.path.abspath(somatic_interactions),
        summary=os.path.abspath(mafsummary), burden_plot=os.path.abspath(burden_plot)
    )

    rmarkdown.render(rmd_script, output_file= os.path.abspath(report_path),
        intermediates_dir = os.path.dirname(os.path.abspath(report_path)),
        output_options=robjects.r.list(self_contained=True), params=parameters
    )

