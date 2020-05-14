from genome_wide_plot import genome_wide_plot
import matplotlib.backends.backend_pdf

# testdata
ideo = "/Users/abramsd/Downloads/cytoBandIdeo.txt"

titan_7 = "/Users/abramsd/work/DATA/QC/titan/Sample007_titan_markers.csv.gz"
roh_7 = "/Users/abramsd/work/DATA/QC/roh/007_roh_ST.txt"
germ_7 = "/Users/abramsd/work/DATA/QC/germline/007/Sample_007_samtools_germline.vcf"
som_7 = "/Users/abramsd/work/DATA/QC/somatic/Sample_007_somatic.csv.gz"
cov_7 = "/Users/abramsd/work/DATA/QC/coverage/merged007_TT"
cov_norm_7 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/7.bins"
remixt_7 = "/Users/abramsd/work/DATA/QC/remixt/results_files_007.h5"
breakpoints_7 = "/Users/abramsd/work/DATA/QC/SVs/Sample_007_filtered_consensus_calls_final.csv.gz"


titan_9 = "/Users/abramsd/work/DATA/QC/titan/Sample009_titan_markers.csv.gz"
roh_9 = "/Users/abramsd/work/DATA/QC/roh/009_roh_ST.txt"
germ_9 = "/Users/abramsd/work/DATA/QC/germline/009/Sample009_samtools_germline.vcf"
som_9 = "/Users/abramsd/work/DATA/QC/somatic/Sample009_somatic.csv.gz"
cov_9 = "/Users/abramsd/work/DATA/QC/coverage/merged009"
cov_norm_9 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_009.tsv"
remixt_9 = "/Users/abramsd/work/DATA/QC/remixt/results_files_009.h5"
breakpoints_9 = "/Users/abramsd/work/DATA/QC/SVs/Sample_009_filtered_consensus_calls_final.csv.gz"

titan_26 = "/Users/abramsd/work/DATA/QC/titan/Sample026_titan_markers.csv.gz"
roh_26 = "/Users/abramsd/work/DATA/QC/roh/026_roh_ST.txt"
germ_26 = "/Users/abramsd/work/DATA/QC/germline/026/Sample_026_samtools_germline.vcf"
som_26 = "/Users/abramsd/work/DATA/QC/somatic/Sample_026_somatic.csv.gz"
cov_26 = "/Users/abramsd/work/DATA/QC/coverage/merged026"
cov_norm_26 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_026.tsv"
remixt_26 = "/Users/abramsd/work/DATA/QC/remixt/results_files_026.h5"
breakpoints_26 = "/Users/abramsd/work/DATA/QC/SVs/Sample_026_filtered_consensus_calls_final.csv.gz"

titan_31 = "/Users/abramsd/work/DATA/QC/titan/Sample031_titan_markers.csv.gz"
roh_31 = "/Users/abramsd/work/DATA/QC/roh/031_roh_ST.txt"
germ_31 = "/Users/abramsd/work/DATA/QC/germline/031/Sample_031_samtools_germline.vcf"
som_31 = "/Users/abramsd/work/DATA/QC/somatic/Sample_031_somatic.csv.gz"
cov_31 = "/Users/abramsd/work/DATA/QC/coverage/merged031"
cov_norm_31 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_031.tsv"
remixt_31 = "/Users/abramsd/work/DATA/QC/remixt/results_Sample_031T"
breakpoints_31 = "/Users/abramsd/work/DATA/QC/SVs/Sample_031_filtered_consensus_calls_final.csv.gz"

titan_36 = "/Users/abramsd/work/DATA/QC/titan/Sample036_titan_markers.csv.gz"
roh_36 = "/Users/abramsd/work/DATA/QC/roh/036_roh_ST.txt"
germ_36 = "/Users/abramsd/work/DATA/QC/germline/036/Sample036_samtools_germline.vcf"
som_36 = "/Users/abramsd/work/DATA/QC/somatic/Sample036_somatic.csv.gz"
cov_36 = "/Users/abramsd/work/DATA/QC/coverage/merged036"
cov_norm_36 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_036.tsv"
remixt_36 = "/Users/abramsd/work/DATA/QC/remixt/results_files_036.h5"
breakpoints_36 = "/Users/abramsd/work/DATA/QC/SVs/Sample_036_filtered_consensus_calls_final.csv.gz"

titan_44 = "/Users/abramsd/work/DATA/QC/titan/Sample044_titan_markers.csv.gz"
roh_44 = "/Users/abramsd/work/DATA/QC/roh/044_roh_ST.txt"
germ_44 = "/Users/abramsd/work/DATA/QC/germline/044/Sample044_samtools_germline.vcf"
som_44 = "/Users/abramsd/work/DATA/QC/somatic/Sample044_somatic.csv.gz"
cov_44 = "/Users/abramsd/work/DATA/QC/coverage/merged044"
cov_norm_44 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_044.tsv"
remixt_44 = "/Users/abramsd/work/DATA/QC/remixt/results_files_044.h5"
breakpoints_44 = "/Users/abramsd/work/DATA/QC/SVs/Sample_044_filtered_consensus_calls_final.csv.gz"

titan_t1 = "/Users/abramsd/work/DATA/QC/titan/Sample_T1_titan_markers.csv.gz"
roh_t1 = "/Users/abramsd/work/DATA/QC/roh/T1_roh_ST.txt"
germ_t1 = "/Users/abramsd/work/DATA/QC/germline/Sample_T1_samtools_germline.vcf"
som_t1 = "/Users/abramsd/work/DATA/QC/somatic/Sample_T1_somatic.csv.gz"
cov_t1 = "/Users/abramsd/work/DATA/QC/coverage/T1T.bins"
cov_norm_t1 = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_T1N.tsv"
remixt_t1 = "/Users/abramsd/work/DATA/QC/remixt/results_T1-T.h5"
breakpoints_t1 = "/Users/abramsd/work/DATA/QC/SVs/Sample_009_filtered_consensus_calls_final.csv.gz"

titan_t2a = "/Users/abramsd/work/DATA/QC/titan/Sample_T2-A_titan_markers.csv.gz"
roh_t2a = "/Users/abramsd/work/DATA/QC/roh/T2A_roh_ST.txt"
germ_t2a = "/Users/abramsd/work/DATA/QC/germline/Sample_T2-A_samtools_germline.vcf"
som_t2a = "/Users/abramsd/work/DATA/QC/somatic/Sample_T2-A_somatic.csv.gz"
cov_t2a = "/Users/abramsd/work/DATA/QC/coverage/T2AT.bins"
cov_norm_t2a = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_T2N.tsv"
remixt_t2a = "/Users/abramsd/work/DATA/QC/remixt/results_T1-T.h5"
breakpoints_t2a = "/Users/abramsd/work/DATA/QC/SVs/Sample_T2-A_filtered_consensus_calls_final.csv.gz"

titan_t2e = "/Users/abramsd/work/DATA/QC/titan/Sample_T2-E_titan_markers.csv.gz"
roh_t2e = "/Users/abramsd/work/DATA/QC/roh/T2E_roh_ST.txt"
germ_t2e = "/Users/abramsd/work/DATA/QC/germline/Sample_T2-E_samtools_germline.vcf"
som_t2e = "/Users/abramsd/work/DATA/QC/somatic/Sample_T2-E_somatic.csv.gz"
cov_t2e = "/Users/abramsd/work/DATA/QC/coverage/T2ET.bins"
cov_norm_t2e = "/Users/abramsd/work/DATA/QC/shipped/normal_coverage/normal_coverage_T2N.tsv"
remixt_t2e = "/Users/abramsd/work/DATA/QC/remixt/results_T2-T-E.h5"
breakpoints_t2e = "/Users/abramsd/work/DATA/QC/SVs/Sample_T2-E_filtered_consensus_calls_final.csv.gz"

titan_002 = "/Users/abramsd/work/DATA/QC/titan/002IO_T_IGO_09443_AQ_7_titan_markers.csv.gz"
roh_002 = "/Users/abramsd/work/DATA/QC/roh/002_roh_ST.txt"
germ_002 = "/Users/abramsd/work/DATA/QC/germline/002IO_T_IGO_09443_AQ_7_samtools_germline.vcf"
som_002 = "/Users/abramsd/work/DATA/QC/somatic/002IO_T_IGO_09443_AQ_7_consensus_somatic.csv.gz"
cov_002 = "/Users/abramsd/work/DATA/QC/coverage/002_T_coverage.tsv"
cov_norm_002 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/002_N_coverage.tsv"
remixt_002 = "/Users/abramsd/work/DATA/QC/remixt/results_002IO_T_IGO_09443_AQ_7"
breakpoints_002 = "/Users/abramsd/work/DATA/QC/SVs/002IO_T_IGO_09443_AQ_7_filtered_consensus_calls_final.csv.gz"

titan_008 = "/Users/abramsd/work/DATA/QC/titan/008LLN_T_IGO_09443_AQ_5_titan_markers.csv.gz"
roh_008 = "/Users/abramsd/work/DATA/QC/roh/008_roh_ST.txt"
germ_008 = "/Users/abramsd/work/DATA/QC/germline/008LLN_T_IGO_09443_AQ_5_samtools_germline.vcf"
som_008 = "/Users/abramsd/work/DATA/QC/somatic/008LLN_T_IGO_09443_AQ_5_consensus_somatic.csv.gz"
cov_008 = "/Users/abramsd/work/DATA/QC/coverage/008_T_coverage.tsv"
cov_norm_008 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/008_N_coverage.tsv"
remixt_008 = "/Users/abramsd/work/DATA/QC/remixt/results_008LLN_T_IGO_09443_AQ_5"
breakpoints_008 = "/Users/abramsd/work/DATA/QC/SVs/008LLN_T_IGO_09443_AQ_5_filtered_consensus_calls_final.csv.gz"

titan_022 = "/Users/abramsd/work/DATA/QC/titan/022LA_T_IGO_09443_AQ_1_titan_markers.csv.gz"
roh_022 = "/Users/abramsd/work/DATA/QC/roh/022_roh_ST.txt"
germ_022 = "/Users/abramsd/work/DATA/QC/germline/022LA_T_IGO_09443_AQ_1_samtools_germline.vcf"
som_022 = "/Users/abramsd/work/DATA/QC/somatic/022LA_T_IGO_09443_AQ_1_consensus_somatic.csv.gz"
cov_022 = "/Users/abramsd/work/DATA/QC/coverage/022_T_coverage.tsv"
cov_norm_022 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/022_N_coverage.tsv"
remixt_022 = "/Users/abramsd/work/DATA/QC/remixt/results_022LA_T_IGO_09443_AQ_1"
breakpoints_022 = "/Users/abramsd/work/DATA/QC/SVs/022LA_T_IGO_09443_AQ_1_filtered_consensus_calls_final.csv.gz"

titan_025 = "/Users/abramsd/work/DATA/QC/titan/025BO_T_IGO_09443_AQ_3_titan_markers.csv.gz"
roh_025 = "/Users/abramsd/work/DATA/QC/roh/025_roh_ST.txt"
germ_025 = "/Users/abramsd/work/DATA/QC/germline/025BO_T_IGO_09443_AQ_3_samtools_germline.vcf"
som_025 = "/Users/abramsd/work/DATA/QC/somatic/0025BO_T_IGO_09443_AQ_3_consensus_somatic.csv.gz"
cov_025 = "/Users/abramsd/work/DATA/QC/coverage/025_T_coverage.tsv"
cov_norm_025 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/025_N_coverage.tsv"
remixt_025 = "/Users/abramsd/work/DATA/QC/remixt/results_025BO_T_IGO_09443_AQ_3"
breakpoints_025 = "/Users/abramsd/work/DATA/QC/SVs/025BO_T_IGO_09443_AQ_3_filtered_consensus_calls_final.csv.gz"

titan_037 = "/Users/abramsd/work/DATA/QC/titan/037RA_T_IGO_09443_AQ_9_titan_markers.csv.gz"
roh_037 = "/Users/abramsd/work/DATA/QC/roh/037_roh_ST.txt"
germ_037 = "/Users/abramsd/work/DATA/QC/germline/037RA_T_IGO_09443_AQ_9_samtools_germline.vcf"
som_037 = "/Users/abramsd/work/DATA/QC/somatic/037RA_T_IGO_09443_AQ_9_consensus_somatic.csv.gz"
cov_037 = "/Users/abramsd/work/DATA/QC/coverage/037_T_coverage.tsv"
cov_norm_037 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/037_N_coverage.tsv"
remixt_037 = "/Users/abramsd/work/DATA/QC/remixt/results_037RA_T_IGO_09443_AQ_9"
breakpoints_037 = "/Users/abramsd/work/DATA/QC/SVs/037RA_T_IGO_09443_AQ_9_filtered_consensus_calls_final.csv.gz"

titan_041 = "/Users/abramsd/work/DATA/QC/titan/041IO_T_IGO_09443_AQ_4_titan_markers.csv.gz"
roh_041 = "/Users/abramsd/work/DATA/QC/roh/041_roh_ST.txt"
germ_041 = "/Users/abramsd/work/DATA/QC/germline/041IO_T_IGO_09443_AQ_4_samtools_germline.vcf"
som_041 = "/Users/abramsd/work/DATA/QC/somatic/041IO_T_IGO_09443_AQ_4_consensus_somatic.csv.gz"
cov_041 = "/Users/abramsd/work/DATA/QC/coverage/041_T_coverage.tsv"
cov_norm_041 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/041_N_NEW_coverage_parsed.tsv"
remixt_041 = "/Users/abramsd/work/DATA/QC/remixt/results_041IO_T_IGO_09443_AQ_4"
breakpoints_041 = "/Users/abramsd/work/DATA/QC/SVs/041IO_T_IGO_09443_AQ_4_filtered_consensus_calls_final.csv.gz"

titan_045 = "/Users/abramsd/work/DATA/QC/titan/045IO_T_IGO_09443_AQ_2_titan_markers.csv.gz"
roh_045 = "/Users/abramsd/work/DATA/QC/roh/045_roh_ST.txt"
germ_045 = "/Users/abramsd/work/DATA/QC/germline/045IO_T_IGO_09443_AQ_2_samtools_germline.vcf"
som_045 = "/Users/abramsd/work/DATA/QC/somatic/045IO_T_IGO_09443_AQ_2_consensus_somatic.csv.gz"
cov_045 = "/Users/abramsd/work/DATA/QC/coverage/045_T_coverage.tsv"
cov_norm_045 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/045_N_coverage.tsv"
remixt_045 = "/Users/abramsd/work/DATA/QC/remixt/results_045IO_T_IGO_09443_AQ_2"
breakpoints_045 = "/Users/abramsd/work/DATA/QC/SVs/045IO_T_IGO_09443_AQ_2_filtered_consensus_calls_final.csv.gz"

titan_049 = "/Users/abramsd/work/DATA/QC/titan/049BO_T_IGO_09443_AQ_8_titan_markers.csv.gz"
roh_049 = "/Users/abramsd/work/DATA/QC/roh/SPECTRUM/049_ST.txt"
germ_049 = "/Users/abramsd/work/DATA/QC/germline/049BO_T_IGO_09443_AQ_8_samtools_germline.vcf"
som_049 = "/Users/abramsd/work/DATA/QC/somatic/049BO_T_IGO_09443_AQ_8_consensus_somatic.csv.gz"
cov_049 = "/Users/abramsd/work/DATA/QC/coverage/049_T_parsed"
cov_norm_049 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/049_N_parsed"
remixt_049 = "/Users/abramsd/work/DATA/QC/remixt/results_049BO_T_IGO_09443_AQ_8"
breakpoints_049 = "/Users/abramsd/work/DATA/QC/SVs/049BO_T_IGO_09443_AQ_8_filtered_consensus_calls_final.csv.gz"

titan_050 = "/Users/abramsd/work/DATA/QC/titan/050IO_T_IGO_09443_AQ_11_titan_markers.csv.gz"
roh_050 = "/Users/abramsd/work/DATA/QC/roh/SPECTRUM/050_ST.txt"
germ_050 = "/Users/abramsd/work/DATA/QC/germline/050IO_T_IGO_09443_AQ_11_samtools_germline.vcf"
som_050 = "/Users/abramsd/work/DATA/QC/somatic/050IO_T_IGO_09443_AQ_11_consensus_somatic.csv.gz"
cov_050 = "/Users/abramsd/work/DATA/QC/coverage/050_T_parsed"
cov_norm_050 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/050_N_parsed"
remixt_050 = "/Users/abramsd/work/DATA/QC/remixt/results_050IO_T_IGO_09443_AQ_11"
breakpoints_050 = "/Users/abramsd/work/DATA/QC/SVs/050IO_T_IGO_09443_AQ_11_filtered_consensus_calls_final.csv.gz"

titan_051 = "/Users/abramsd/work/DATA/QC/titan/051IO_T_IGO_09443_AQ_12_titan_markers.csv.gz"
roh_051 = "/Users/abramsd/work/DATA/QC/roh/SPECTRUM/051_ST.txt"
germ_051 = "/Users/abramsd/work/DATA/QC/germline/051IO_T_IGO_09443_AQ_12_samtools_germline.vcf"
som_051 = "/Users/abramsd/work/DATA/QC/somatic/051IO_T_IGO_09443_AQ_12_consensus_somatic.csv.gz"
cov_051 = "/Users/abramsd/work/DATA/QC/coverage/051_T_parsed"
cov_norm_051 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/051_N_parsed"
remixt_051 = "/Users/abramsd/work/DATA/QC/remixt/results_051IO_T_IGO_09443_AQ_12"
breakpoints_051 = "/Users/abramsd/work/DATA/QC/SVs/051IO_T_IGO_09443_AQ_12_filtered_consensus_calls_final.csv.gz"

titan_053 = "/Users/abramsd/work/DATA/QC/titan/053IO_T_IGO_09443_AQ_14_titan_markers.csv.gz"
roh_053 = "/Users/abramsd/work/DATA/QC/roh/SPECTRUM/053_ST.txt"
germ_053 = "/Users/abramsd/work/DATA/QC/germline/053IO_T_IGO_09443_AQ_14_samtools_germline.vcf"
som_053 = "/Users/abramsd/work/DATA/QC/somatic/053IO_T_IGO_09443_AQ_14_consensus_somatic.csv.gz"
cov_053 = "/Users/abramsd/work/DATA/QC/coverage/053_T_parsed"
cov_norm_053 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/053_N_parsed_SWAPPED"
remixt_053 = "/Users/abramsd/work/DATA/QC/remixt/results_053IO_T_IGO_09443_AQ_14"
breakpoints_053 = "/Users/abramsd/work/DATA/QC/SVs/053IO_T_IGO_09443_AQ_14_filtered_consensus_calls_final.csv.gz"

titan_054 = "/Users/abramsd/work/DATA/QC/titan/054IO_T_IGO_09443_AQ_13_titan_markers.csv.gz"
roh_054 = "/Users/abramsd/work/DATA/QC/roh/SPECTRUM/054_ST.txt"
germ_054 = "/Users/abramsd/work/DATA/QC/germline/054IO_T_IGO_09443_AQ_13_samtools_germline.vcf"
som_054 = "/Users/abramsd/work/DATA/QC/somatic/054IO_T_IGO_09443_AQ_13_consensus_somatic.csv.gz"
cov_054 = "/Users/abramsd/work/DATA/QC/coverage/054_T_parsed"
cov_norm_054 = "/Users/abramsd/work/DATA/QC/normal_bam_coverage/054_N_parsed_SWAPPED"
remixt_054 = "/Users/abramsd/work/DATA/QC/remixt/results_054IO_T_IGO_09443_AQ_13"
breakpoints_054 = "/Users/abramsd/work/DATA/QC/SVs/054IO_T_IGO_09443_AQ_13_filtered_consensus_calls_final.csv.gz"


# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_007.pdf")
# genome_wide_plot(remixt_7, "007", titan_7, roh_7, germ_7, som_7, cov_7, cov_norm_7, breakpoints_7, ideo,["1", "4"],
#                  "/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_007.pdf")
# # pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_009.pdf")
# genome_wide_plot(remixt_9, "009", titan_9, roh_9, germ_9, som_9, cov_9, cov_norm_9, breakpoints_9, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_026.pdf")
# genome_wide_plot(remixt_26, "026", titan_26, roh_26, germ_26, som_26, cov_26, cov_norm_26, breakpoints_26, ideo, pdf)
# pdf.close()

# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_031.pdf")
# # genome_wide_plot(remixt_31, "031", titan_31, roh_31, germ_31, som_31, cov_31, cov_norm_31, breakpoints_31, ideo, pdf)
# # pdf.close()

# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_036.pdf")
# genome_wide_plot(remixt_36, "036", titan_36, roh_36, germ_36, som_36, cov_36, cov_norm_36, breakpoints_36, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_044.pdf")
# genome_wide_plot(remixt_44, "044", titan_44, roh_44, germ_44, som_44, cov_44, cov_norm_44, breakpoints_44, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_T1.pdf")
# genome_wide_plot(remixt_t1, "T1", titan_t1, roh_t1, germ_t1, som_t1, cov_t1, cov_norm_t1, breakpoints_t1, ideo, pdf)
# pdf.close()
# #
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_T2-A.pdf")
# genome_wide_plot(remixt_t2a, "T2-A", titan_t2a, roh_t2a, germ_t2a, som_t2a, cov_t2a, cov_norm_t2a, breakpoints_t2a, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_T2-E.pdf")
# genome_wide_plot(remixt_t2e, "T2-E", titan_t2e, roh_t2e, germ_t2e, som_t2e, cov_t2e, cov_norm_t2e, breakpoints_t2e, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_002.pdf")
# genome_wide_plot(remixt_002, "002", titan_002, roh_002, germ_002, som_002, cov_002, cov_norm_002, breakpoints_002, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_008.pdf")
# genome_wide_plot(remixt_008, "008", titan_008, roh_008, germ_008, som_008, cov_008, cov_norm_008, breakpoints_008, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_022.pdf")
# genome_wide_plot(remixt_022, "022", titan_022, roh_022, germ_022, som_022, cov_022, cov_norm_022, breakpoints_022, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_025.pdf")
# genome_wide_plot(remixt_025, "025", titan_025, roh_025, germ_025, som_025, cov_025, cov_norm_025, breakpoints_025, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_titan_037.pdf")
# genome_wide_plot(remixt_037, "037", titan_037, roh_037, germ_037, som_037,
#                  cov_037, cov_norm_037, breakpoints_037, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_041.pdf")
# genome_wide_plot(remixt_041, "041", titan_041, roh_041, germ_041, som_041, cov_041, cov_norm_041, breakpoints_041, ideo, pdf)
# pdf.close()
#
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_045.pdf")
# genome_wide_plot(remixt_045, "045", titan_045, roh_045, germ_045, som_045, cov_045, cov_norm_045, breakpoints_045, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_049.pdf")
# genome_wide_plot(remixt_049, "049", titan_049, roh_049, germ_049, som_049, cov_049, cov_norm_049, breakpoints_049, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_050.pdf")
# genome_wide_plot(remixt_050, "050", titan_050, roh_050, germ_050, som_050, cov_050, cov_norm_050, breakpoints_050, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_051.pdf")
# genome_wide_plot(remixt_051, "051", titan_051, roh_051, germ_051, som_051, cov_051, cov_norm_051, breakpoints_051, ideo, pdf)
# pdf.close()
#
# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_053.pdf")
# genome_wide_plot(remixt_053, "053", titan_053, roh_053, germ_053, som_053, cov_053, cov_norm_053, breakpoints_053, ideo, pdf)
# pdf.close()

# pdf = matplotlib.backends.backend_pdf.PdfPages("/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_054.pdf")
genome_wide_plot(remixt_054, "054", titan_054, roh_054, germ_054, som_054, cov_054, cov_norm_054, breakpoints_054, ideo, ["8"], "/Users/abramsd/work/OUTPUTS/QC/genome_wide_remixt/genome_wide_remixt_054.pdf")
# pdf.close()

# from wgs_qc_utils.reader import read_remixt
# from wgs_qc_utils.plotter import remixt_plotting
# import matplotlib.pyplot as plt
# remixt = read_remixt.read(remixt_054, "054")
# remixt = remixt_plotting.prepare_at_chrom(remixt, "1")
# fig, ax = plt.subplots()
# ax = remixt_plotting.plot(remixt, ax, remixt.start.max()/1000000)
# plt.show()
