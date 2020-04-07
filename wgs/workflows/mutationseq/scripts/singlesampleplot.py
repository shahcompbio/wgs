import matplotlib
import pandas as pd
import pysam

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse

from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np


class PlotSingleSample(object):
    def __init__(self, variants, thousand_genomes, dbsnp, concordance_ref, outpdf, threshold):
        self.variants = variants
        self.thousand_genomes = thousand_genomes
        self.dbsnp = dbsnp
        self.concordance_ref = concordance_ref
        self.threshold = threshold

        self.outpdf = outpdf

        self.total_calls_in_ref = 0

    def load_concordance_reference(self, tablename):
        with pd.HDFStore(self.concordance_ref, 'r') as indata:
            return indata[tablename][0].to_dict()

    def load_reference_data(self, chromosome):
        positions = set()

        for filename in (self.thousand_genomes, self.dbsnp):
            with pysam.Tabixfile(filename) as ref_vcf_data:
                for record in ref_vcf_data.fetch(chromosome):
                    record = record.replace('chr', '').split()

                    k = record[0] + ',' + record[1]
                    positions.add(k)

        return positions

    def get_variant_prob_counts(self):

        all_pos = defaultdict(int)
        matched_pos = defaultdict(int)

        with pysam.Tabixfile(self.variants) as variant_data:
            chromosomes = variant_data.contigs
            for chromosome in chromosomes:
                ref_pos = self.load_reference_data(chromosome)

                self.total_calls_in_ref += len(ref_pos)

                for record in variant_data.fetch(chromosome):
                    record = record.replace('chr', '').split()

                    k = record[0] + ',' + record[1]
                    pr = float(record[7].split(';')[0].split('=')[1])

                    if float(pr) >= self.threshold:
                        all_pos[pr] += 1
                        if k in ref_pos:
                            matched_pos[pr] += 1

        return matched_pos, all_pos

    def convert_to_cumulative(self, counts):
        sorted_keys = sorted(counts.keys(), reverse=True)
        values = np.cumsum([counts[key] for key in sorted_keys])
        cum_counts = {k: v for k, v in zip(sorted_keys, values)}
        return cum_counts

    def get_percentages(self, counts_input, counts_ref, cumulative=False):
        if cumulative:
            counts_input = self.convert_to_cumulative(counts_input)
            counts_ref = self.convert_to_cumulative(counts_ref)

        keys = set(counts_input.keys()).intersection(set(counts_ref.keys()))
        keys = sorted(keys)
        data = {}
        for k in keys:
            data[k] = float(counts_ref[k]) / counts_input[k] * 100
        return data

    def concordance_plot(self, counts_all, counts_ref, subplot_axes, plot):
        counts_all = self.convert_to_cumulative(counts_all)
        counts_ref = self.convert_to_cumulative(counts_ref)

        percentages = self.get_percentages(counts_all, counts_ref)
        ref_percentages = self.load_concordance_reference('concordance_percentage')

        ax = plot.add_subplot(subplot_axes[0], subplot_axes[1], subplot_axes[2])
        ax.set_title('concordance', fontsize=8)
        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("percentage", fontsize=8)

        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)

        keys = sorted(percentages.keys())
        ax.plot(keys, [percentages[k] for k in keys], linewidth=1.5, label='concordance curve')

        keys = sorted(ref_percentages.keys())
        ax.plot(keys, [ref_percentages[k] for k in keys], linewidth=1.5, label='concordance reference')

        plt.ylim(0, 100)
        ax.legend(loc='lower right', fontsize=6, handlelength=5)

    def counts_plot(self, counts_all, counts_matched, subplot_axes, plot, loc='upper right', cumulative=False):

        if cumulative:
            counts_all = self.convert_to_cumulative(counts_all)
            counts_matched = self.convert_to_cumulative(counts_matched)
            counts_all_ref = self.load_concordance_reference('cumulative_all_counts')
            counts_matched_ref = self.load_concordance_reference('cumulative_matched_counts')
            title = 'Number of snvs for all positions (solid) and \n positions found (Dashed)'
        else:
            counts_all_ref = self.load_concordance_reference('absolute_all_counts')
            counts_matched_ref = self.load_concordance_reference('absolute_matched_counts')
            title = 'Density of all (solid) and matched (Dashed) positions'

        ax = plot.add_subplot(subplot_axes[0], subplot_axes[1], subplot_axes[2])

        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)

        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)

        keys = sorted(counts_all.keys())
        ax.plot(keys, [counts_all[k] for k in keys], linestyle='solid', label='all positions')

        keys = sorted(counts_matched.keys())
        ax.plot(keys, [counts_matched[k] for k in keys], linestyle='dashed', label='matched positions')

        keys = sorted(counts_all_ref.keys())
        ax.plot(keys, [counts_all_ref[k] for k in keys], linestyle='solid', label='ref all positions')

        keys = sorted(counts_matched_ref.keys())
        ax.plot(keys, [counts_matched_ref[k] for k in keys], linestyle='dashed', label='ref matched positions')

        ax.legend(loc=loc, fontsize=4, handlelength=5)

        ax.set_title(title, fontsize=8)

    def sensitivity_plot(self, counts_matched, subplot_axes, plot):

        counts_matched = self.convert_to_cumulative(counts_matched)
        sensitivity_ref = self.load_concordance_reference('sensitivity')

        ax = plot.add_subplot(subplot_axes[0], subplot_axes[1], subplot_axes[2])
        ax.set_title('sensitivity', fontsize=8)
        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("Sensitivity(%)", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)

        keys_ref = sorted(set(counts_matched.keys()))
        sensitivity = []
        for k in keys_ref:
            sensitivity.append(float(counts_matched.get(k)) / self.total_calls_in_ref)

        ax.plot(list(keys_ref), sensitivity, linewidth=1.5, label='sensitivity curve')
        ax.legend(loc='lower left', fontsize=6, handlelength=5)

        keys_ref = sorted(sensitivity_ref.keys())
        sensitivity = [sensitivity_ref[k] for k in keys_ref]
        ax.plot(list(keys_ref), sensitivity, linewidth=1.5, label='reference curve')
        ax.legend(loc='lower left', fontsize=6, handlelength=5)

    def get_top3_counts(self, counts_all, counts_ref):
        counts_all = self.convert_to_cumulative(counts_all)
        counts_ref = self.convert_to_cumulative(counts_ref)

        value_closest_3m = min(counts_all.values(), key=lambda x: abs(x - 3e6))
        pr = [pr for pr, count in counts_all.items() if count == value_closest_3m][0]
        concordance_3m = self.get_percentages(counts_all, counts_ref)[pr]

        value_closest_35m = min(counts_all.values(), key=lambda x: abs(x - 3.5e6))
        pr = [pr for pr, count in counts_all.items() if count == value_closest_35m][0]
        concordance_35m = self.get_percentages(counts_all, counts_ref)[pr]

        return concordance_3m, value_closest_3m, concordance_35m, value_closest_35m

    def plot_top_con(self, counts_all, counts_ref, subplot_axes, plot):

        top3_con, top3_val, top35_con, top35_val = self.get_top3_counts(counts_all, counts_ref)

        ref_top_concordance = self.load_concordance_reference('top_concordance')
        ref_top_concordance_3m = ref_top_concordance['concordance_top3_mil']
        ref_top_concordance_35m = ref_top_concordance['concordance_top35_mil']

        label_3 = 'top ' + str(top3_val) + ' positions'
        label_35 = 'top ' + str(top35_val) + ' positions'

        ax = plot.add_subplot(subplot_axes[0], subplot_axes[1], subplot_axes[2])

        ax.bar(1, ref_top_concordance_3m, label='reference 3e6 positions')
        ax.bar(3, ref_top_concordance_35m, label='reference 35e5 positions')
        ax.bar(0, top3_con, label=label_3)
        ax.bar(2, top35_con, label=label_35)

        ax.set_title('Concordance for top 3/3.5 million positions', fontsize=8)

        ax.set_ylabel("Concordance", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)

        ax.xaxis.set_tick_params(labelsize=0)
        ax.legend(loc='lower right',
                  fontsize=6)

    def main(self):
        matched_counts, all_counts = self.get_variant_prob_counts()

        pdfout = PdfPages(self.outpdf)
        plot = plt.figure(figsize=(10, 10))
        plot.suptitle('Mutationseq (Single Sample) portraits', fontsize=9)

        self.concordance_plot(all_counts, matched_counts, (3, 2, 1), plot)
        self.counts_plot(all_counts, matched_counts, (3, 2, 2), plot, cumulative=True)
        self.counts_plot(all_counts, matched_counts, (3, 2, 4), plot, cumulative=False)

        self.sensitivity_plot(matched_counts, (3, 2, 3), plot)

        self.plot_top_con(all_counts, matched_counts, (3, 2, 5), plot)

        plt.tight_layout()

        pdfout.savefig(plot)
        pdfout.close()


def parse_args():
    parser = argparse.ArgumentParser(prog='generate feature boxplots and concordance plots',
                                     description='''generates boxplots for feature
                                                    distributions and concordance plots''')

    parser.add_argument('-v', '--variant_file',
                        required=True,
                        help='Path to the variant file (/path/to/vcf)')

    parser.add_argument('--output',
                        required=True,
                        help='''path to the folder with prefix for the
                                outputs (/path/to/output_folder/name)''')

    parser.add_argument('--dbsnp',
                        help='path to dbsnp reference file')

    parser.add_argument('--thousand_gen',
                        help='path to thousand genomes reference file')

    parser.add_argument('--threshold',
                        default=0.5,
                        type=float,
                        help=''' threshold''')

    parser.add_argument("-d", "--ref_data",
                        default=None,
                        help='''Adds reference information to plots''')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    ssplt = PlotSingleSample(
        args.variant_file,
        args.thousand_gen,
        args.dbsnp,
        args.ref_data,
        args.output,
        args.threshold)
    ssplt.main()
