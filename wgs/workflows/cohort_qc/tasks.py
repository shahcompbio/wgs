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
    '''
    make list of genes to display in plots
    Parameters
    ----------
    cna_table : (unused) copynumber data
    genelist : custom list of genes to display
    base_genes : default = default set of genes to display
    Returns
    -------
    '''
    # TODO: derive gene list from cna data emperically 
    if base_genes == 'default':
        genes = ["PPM1D", "TP53", "BRCA1", "BRCA2", "MECOM", "RB1", "PTEN", "PALB2",
                 "ERBB2", "CDK12", "PIK3CA", "KRAS", "CCNE1", "MYC"]
    genes = pd.DataFrame({"gene": list(genes)})
    genes.to_csv(genelist)


def merge_segmental_cn(segmental_cn, concats):
    '''
    merge segmental copynumber files into concats
    Parameters
    ----------
    segmental_cn : dictionary of filepaths
    concats : merged filepath
    Returns
    -------
    '''
    files = [pd.read_csv(f, sep="\t") for f in list(segmental_cn.values())]
    segmental_cn_combined = pd.concat(files)
    segmental_cn_combined.to_csv(concats, sep="\t", index=False)


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

    segment_cn(segmental_cn, aggregated_cn_data, stats["ploidy"])


def segment_cn(filename, aggregated_cn_data, ploidy,  
        cn_col="copy", length_col="length"
    ):
    '''
    segment copynumber
    Parameters
    ----------
    filename : output path
    aggregated_cn_data : copynumber data
    Returns
    -------
    '''
    aggregated_cn_data['ploidy'] = ploidy
    aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data[cn_col] 
        / aggregated_cn_data['ploidy']
    )
    aggregated_cn_data['num.mark'] = (aggregated_cn_data[length_col] / 500000).astype(int)
    aggregated_cn_data = aggregated_cn_data.rename(
        columns={'sample': 'ID', 'chromosome': 'chrom', 
            'start': 'loc.start', 'end': 'loc.end'
        }
    )
    aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 
        'loc.end', 'num.mark', 'seg.mean'
    ]]
    aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
    aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
    aggregated_cn_data = transformations._correct_seg_bin_ends(aggregated_cn_data)
    aggregated_cn_data.to_csv(filename, index=None, sep='\t')


def make_cbio_cna_table(cn_change_filename, cbio_table):
    '''
    make copynumber cbioportal input
    Parameters
    ----------
    cn_change_filename : input cn
    cbio_table : output path
    Returns
    -------
    '''
    gistic_data = pd.read_csv(cn_change_filename, sep="\t", 
        usecols=["gene_name", "sample", "gistic_value"])
    gistic_data = gistic_data.rename(columns={"gene_name": "Hugo_Symbol"})
    gistic_data = gistic_data.drop_duplicates()
    gistic_data = gistic_data.astype({"gistic_value": "Int64"})

    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)
    gistic_matrix.to_csv(cbio_table, sep="\t", index=False, na_rep="NA")


def make_maftools_cna_table(cn_change_filename, maftools_table):
    '''
    format cn data for maftools
    Parameters
    ----------
    cn_change_filename : input cn
    maftools_table : output path
    Returns
    -------
    '''
    cn_change = pd.read_csv(cn_change_filename, sep="\t",
                            usecols=["gene_name", "sample", "is_hdel", "is_loh", "is_hlamp"])

    cn_change = cn_change[cn_change['is_hdel'] | cn_change['is_loh'] | cn_change['is_hlamp']]

    cn_change["is_loh"] = cn_change.is_loh.replace({False: "", True: "loh"})
    cn_change["is_hdel"] = cn_change.is_hdel.replace({False: "", True: "hdel"})
    cn_change["is_hlamp"] = cn_change.is_hlamp.replace({False: "", True: "hlamp"})
    cn_change["CN"] = cn_change[["is_hdel", "is_loh", "is_hlamp"]].apply(
        lambda row: "".join(map(str, row)) ,axis=1
    )

    cn_change = cn_change.rename(columns={"gene_name": "Gene", "cn_type": 
        "CN", "sample": "Sample_Name"})
    cn_change = cn_change[["Gene", "Sample_Name", "CN"]]
    cn_change.to_csv(maftools_table, index=False, sep="\t")


def merge_cna_tables(tables, output):
    '''
    merge copynumber tables
    Parameters
    ----------
    tables : dictionary of paths
    output : output path
    Returns
    -------
    '''
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


def classify_remixt(sample_label, remixt, gtf, cn_change):
    '''
    run classifycopynumber of remixt data
    Parameters
    ----------
    remixt : remixt postprocessing cn file
    sample_label : sample label for remixt
    gtf: gtf file to use in classifycopynumber
    cn_change: cn change file
    Returns
    -------
    '''
    cmd = [
        "classifycopynumber", gtf, cn_change, "--remixt_parsed_csv", 
        remixt, '--sample_ids', sample_label]
    pypeliner.commandline.execute(*cmd)


def _write_maf(m, label, merged_maf, write_header):
    '''
    write a maf to file appending
    Parameters
    ----------
    m : maf file
    label : file label
    merged_maf: merged output
    write_header: T/F
    Returns
    -------
    '''
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

        chunk.to_csv(merged_maf, sep="\t", index=False, header=write_header, 
            mode='a', na_rep="")
        write_header=False


def merge_mafs(germline, somatic_mafs, merged_maf):
    '''
    merge a set of mafs into one file
    Parameters
    ----------
    germline : germline maf path dictionary 
    somatic_mafs : somatic maf path dictionary
    merged_maf: merged output
    Returns
    -------
    '''
    write_header = True

    for label, m in somatic_mafs.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header = False
    for label, m in germline.items():
        _write_maf(m, label, merged_maf, write_header)
        write_header = False


def annotate_maf_with_oncokb(maf, api_key, tmpspace, annotated_maf):
    '''
    annotate maf with oncokb 
    Parameters
    ----------
    maf :  maf path to annotate 
    somatic_mafs : somatic maf path dictionary
    merged_maf: merged output
    Returns
    -------
    '''
    helpers.makedirs(tmpspace)
    ma.annotate(maf, annotated_maf, api_key)


def filter_maf(annotated_maf, filtered_maf, write_header=True):
    '''
    filter maf based on annotated
    Parameters
    ----------
    annotated_maf :  maf annotated with oncokb
    filtered_maf : output path
    write_header: T/F
    Returns
    -------
    '''
    oncogenic_annotations = ["Oncogenic", "Likely Oncogenic", "Predicted Oncogenic"]
    maf = pd.read_csv(annotated_maf, sep="\t", dtype='str', chunksize=10e6)
    for chunk in maf:
        chunk = chunk[chunk.ONCOGENIC.isin(oncogenic_annotations)]
        chunk.to_csv(filtered_maf, sep="\t", index=False, header=write_header, mode='a')

        write_header = False


def annotate_germline_somatic(filtered_maf, annotated_maf, is_germline):
    '''
    annotate maf entries with germline/somatic info
    Parameters
    ----------
    filtered_maf :  filtered mf
    annotated_maf : output annotated maf
    is_germline: T/F
    Returns
    -------
    '''
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
    maf["is_germline"] = [is_germline] * len(maf)
    maf.to_csv(annotated_maf, sep="\t", index=False)


def label_germline_somatic(row):
    '''
    label a maf row germline/somatic
    Parameters
    ----------
    row :  pandas row
    Returns
    -------
    '''
    if row.is_germline == True:
        return row.Variant_Classification + "_" + "germline"
    return row.Variant_Classification + "_" + "somatic"


def prepare_maf_for_maftools(cohort_label, filtered_maf, prepared_maf, 
        non_synonymous_labels, vcNames):
    '''
    filter on non synonymous labels
    add germline/somatic annotate to variant classification
    add group/patient labels
    write out vcNames
    '''
    maf = pd.read_csv(filtered_maf, sep="\t", dtype='str')
    maf = maf[maf.Variant_Classification.isin(non_synonymous_labels)]
    maf["Variant_Classification"] = maf.apply(lambda row: label_germline_somatic(row), axis=1)
    nonsynclasses = pd.DataFrame(
        {"Variant_Classification": 
            maf.Variant_Classification.unique().tolist()
        }
    )

    nonsynclasses.to_csv(vcNames, index=False)
    maf.to_csv(prepared_maf, sep="\t", index=False)


def plot_mutation_burden(maf, burden_plot_path):
    '''
    plot the mutation burden for samples in the maf
    Parameters
    ----------
    maf :  annotated, filtered maf
    burden_plot_path : output path for plot
    Returns
    -------
    '''
    maf = pd.read_csv(maf, dtype='str', sep="\t", 
        usecols=["Tumor_Sample_Barcode"]).drop_duplicates()
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
        mafsummary, vcNames, genelist
):
    '''
    make maftools plots for maf 
    Parameters
    ----------
    maf :  annotated, filtered maf
    cntable : copynumber data
    oncoplot_path: path of oncoplot
    somatic_interactions : path of somatic plot
    mafsummary : path of maf summary plot
    vcNames : set of alteration types in maf
    geneList : list of genes to show in plots
    Returns
    -------
    '''
    maftools_plots_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'scripts','maftools_plots.R'
    )
    plots_cmd = [
        "Rscript", maftools_plots_script, maf, vcNames, cntable,
        oncoplot_path, somatic_interactions, mafsummary, genelist
    ]

    pypeliner.commandline.execute(*plots_cmd)



def make_report(cohort_label, oncoplot, somatic_interactions, 
    mafsummary, burden_plot, report_path
):
    '''
    make report from plots
    Parameters
    ----------
    cohort_label :  label of cohort
    oncoplot_path: path of oncoplot
    somatic_interactions : path of somatic plot
    mafsummary : path of maf summary plot
    burden_plot : path of burden plot
    report_path : path to report
    Returns
    -------
    '''
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

