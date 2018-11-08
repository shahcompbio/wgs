'''
Created on Apr 2, 2014

@author: dgrewal
@modified: 5 Feb 2015 by jrosner (self.args.reference -> dbsnp and thousand_gen)
'''
import logging
import os
import pysam
import subprocess
import re
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class concordance_plots(object):
    def __init__(self,args):
        self.args = args
        self.colors = ['b', 'm', 'g', 'y', 'k', 'c']
        self.plot = None
        self.label = None
        # matplotlib.rcParams['axes.color_cycle'] = ['b', 'r', 'g', 'y', 'k', 'c', 'm']
        self.outfile = open(self.args.data,'w')

        self.len_ref = 0
        self.cum_ref_counts = defaultdict(int)
        self.cum_all_counts = defaultdict(int)
        self.abs_ref_counts = defaultdict(int)
        self.abs_all_counts = defaultdict(int)
        self.sensitivity = defaultdict(int)

        ref_labs = []
        if self.args.dbsnp:
            ref_labs = ['dbsnp']

        if self.args.thousand_gen:
            ref_labs.append('1000Gen')

        self.ref_labels = '+'.join(ref_labs)
        self.references = [x for x in [self.args.dbsnp, self.args.thousand_gen] if x is not None]

        self.variant_label = self.args.variant_label
        self.chromosomes = map(str,range(1,23))
        self.chromosomes.append('X')
        self.chromosomes.append('Y')
        self.delfiles = []
        self.__check_inputs()
        self.con_35_ref = True

        self.ref_con_perc = None
        self.ref_cum_all_counts = None
        self.ref_cum_mat_counts = None
        self.ref_con_3mil = None
        self.ref_con_35mil = None
        self.ref_abs_all_counts = None
        self.ref_abs_mat_counts = None
        self.ref_sens = None
        self.museq_percentages = None


    def __tabix_file(self,filename):

        sort_cmd = "cat {} | vcf-sort > {}.sorted.vcf"
        gz_cmd = 'bgzip -f -c '+ filename +' > '+ filename+'.gz'
        tbx_cmd = 'tabix '+ filename+'.gz'



        logging.info('bgzip and tabix '+filename)
        gz_p = subprocess.Popen(gz_cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell = True)
        out = gz_p.communicate()
        if out != ('',''):
            logging.error('error gzipping file:' + filename+ 'error:'+str(out))
            raise Exception('error gzipping file:' + filename+ 'error:'+str(out))

        tbx_p = subprocess.Popen(tbx_cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell = True)
        out = tbx_p.communicate()

        #list the file for deletion
        self.delfiles.append(filename+'.gz')
        self.delfiles.append(filename+'.gz.tbi')

        if out != ('',''):
            logging.error('error tabix\'ing file:' + filename+ 'error:'+str(out))
            raise Exception('error tabix\'ing file:' + filename+ 'error:'+str(out))

    def __is_vcf(self,filename):
        file_str = open(filename)
        for line in file_str:
            if re.match(r'^\w+$|#|=', line.strip().split()[0]):
                return True
            else:
                return False

    def __check_inputs(self):
        new_ref = []
        for ref in self.references:
            if not self.__is_vcf(ref):
                new_ref.append(ref)
            else:
                self.__tabix_file(ref)
                new_ref.append(ref+'.gz')

        self.references = new_ref

        #check and index variant file
        if self.__is_vcf(self.args.variant_file):
            self.__tabix_file(self.args.variant_file)
            self.args.variant_file = self.args.variant_file+'.gz'


    def __parse_ref(self,ref_set,chromosome,filename):
        try:
            tabix_file = pysam.Tabixfile(filename)
            for record in tabix_file.fetch(chromosome):
                record = record.replace('chr','').split()

                k = record[0]+','+record[1]
                ref_set.add(k)
        except:
            #but continue execution
            logging.error('chr '+str(chromosome)+' missing in: '+filename)


    def __get_ref_positions(self,chromosome):
        logging.info('Generating the reference set for chromosome: '+chromosome)
        ref_set = set()

        for ref in self.references:
            filename = ref.strip()
            self.__parse_ref(ref_set,chromosome,filename)

        if len(self.references) == 0:
            logging.error('At least one of references should be provided')
            raise Exception('Please input reference files (dbsnp and/or 1000 genomes)')

        return ref_set

    def __read_museq_vcf(self,chromosome):
        ref_set = self.__get_ref_positions(chromosome)
        self.len_ref += len(ref_set)

        logging.info('Filtering the museq vcf file, chromosome:'+chromosome)

        matched_pos = set()
        all_pos = set()

        try:
            file_name = self.args.variant_file.strip()
            tabix_file = pysam.Tabixfile(file_name)
            for record in tabix_file.fetch(chromosome):
                record = record.replace('chr','').split()

                k = record[0]+','+record[1]
                pr = record[7].split(';')[0].split('=')[1]

                if float(pr)>=self.args.threshold:
                    all_pos.add((k,pr))
                    if k in ref_set:
                        matched_pos.add((k,pr))
            ref_set = []
        except:
            #but continue execution
            logging.error('chr '+str(chromosome)+' missing in: '+file_name)

        return matched_pos,all_pos


    #get the counts for each probability
    def __get_prob_set(self,input_set,threshold):
        probabilities = []
        if threshold == self.args.threshold:
            probabilities = [pr[1] for pr in input_set]
        else:
            for val in input_set:
                pr = val[1]
                if pr >= threshold:
                    probabilities.append(pr)
        return probabilities

    def __get_counts(self,probabilities,cumulative):
        counts = defaultdict(int)

        for i in probabilities:
            counts[i] +=1

        #positions with p>=given prob
        if cumulative:
            keys = sorted(counts.keys(), reverse=True)
            x = [sum([counts[keys[j]] for j in range(i)]) for i in range(len(keys))]
            for i in range(len(keys)):
                k = keys[i]
                counts[k] += x[i]
        return counts

    def __get_percentages(self,counts_input,counts_ref):
        keys = set(counts_input.keys()).intersection(set(counts_ref.keys()))
        keys = sorted(keys)
        data = {}
        for k in keys:
            data[k] = float(counts_ref[k])/counts_input[k] * 100
        return data

    def __merge_counts_dict(self,cum_counts_all,cum_counts_ref,abs_counts_all,abs_counts_ref):

        for key,value in cum_counts_all.iteritems():
            if self.cum_all_counts.get(key) == None:
                self.cum_all_counts[key]=0
            self.cum_all_counts[key] = self.cum_all_counts.get(key)+value

        for key,value in cum_counts_ref.iteritems():
            if self.cum_ref_counts.get(key) == None:
                self.cum_ref_counts[key]=0
            self.cum_ref_counts[key] = self.cum_ref_counts.get(key)+value

        for key,value in abs_counts_all.iteritems():
            if self.abs_all_counts.get(key) == None:
                self.abs_all_counts[key]=0
            self.abs_all_counts[key] = self.abs_all_counts.get(key)+value

        for key,value in abs_counts_ref.iteritems():
            if self.abs_ref_counts.get(key) == None:
                self.abs_ref_counts[key]=0
            self.abs_ref_counts[key] = self.abs_ref_counts.get(key)+value


    def __get_data(self):
        for chrom in self.chromosomes:
            matched_pos,all_pos = self.__read_museq_vcf(chrom)
            input_prob = self.__get_prob_set(all_pos,self.args.threshold)
            ref_prob = self.__get_prob_set(matched_pos,self.args.threshold)
            cum_counts_all = self.__get_counts(input_prob,True)
            cum_counts_ref = self.__get_counts(ref_prob,True)
            abs_counts_all = self.__get_counts(input_prob,False)
            abs_counts_ref = self.__get_counts(ref_prob,False)
            self.__merge_counts_dict(cum_counts_all, cum_counts_ref,abs_counts_all,abs_counts_ref)

    def __read_ref_data(self):
        file_str = open(self.args.ref_data)
        for line in file_str:
            line = line.strip().split()
            if line[0] != '#':
                logging.error('Error in the reference data file')
                return None

            if line[1] == 'concordance_percentage':
                self.ref_con_perc = eval(line[2])

            elif line[1] =='cumulative_all_counts':
                self.ref_cum_all_counts = eval(line[2])

            elif line[1] =='cumulative_matched_counts':
                self.ref_cum_mat_counts = eval(line[2])

            elif line[1] =='concordance_top3_mil':
                self.ref_con_3mil = float(line[2])

            elif line[1] =='concordance_top35_mil':
                self.ref_con_35mil = float(line[2])

            elif line[1] =='absolute_all_counts':
                self.ref_abs_all_counts = eval(line[2])

            elif line[1] =='absolute_matched_counts':
                self.ref_abs_mat_counts = eval(line[2])

            elif line[1] =='sensitivity':
                self.ref_sens = eval(line[2])

    def __plot_ref_data(self):
        self.__concordance_plot(self.ref_con_perc,(3,2,1), 'reference curve')
        self.__counts_plot(self.ref_cum_all_counts,(3,2,2),'reference all positions',False)
        self.__counts_plot(self.ref_cum_mat_counts,(3,2,2),'reference matched positions',True)
        self.__sensitivity_plot((3,2,3),'reference curve',sensitivity = self.ref_sens)

        self.__counts_plot(self.ref_abs_all_counts,(3,2,4),'reference all positions',False,'upper left')
        self.__counts_plot(self.ref_abs_mat_counts,(3,2,4),'reference matched positions',True,'upper left')

        if self.con_35_ref:
            self.__plot_top_con(self.ref_con_3mil,'3000000',self.ref_con_35mil,'3500000',(3,2,5),ref = True)


    def __generate_plots(self):
        logging.info('Calculating concordance and frequencies... ')

        self.__get_data()
        self.museq_percentages = self.__get_percentages(self.cum_all_counts,self.cum_ref_counts)

        self.__concordance_plot(self.museq_percentages,(3,2,1),'concordance curve')
        self.__counts_plot(self.cum_all_counts,(3,2,2),'all positions',False)
        self.__counts_plot(self.cum_ref_counts,(3,2,2),'matched positions',True)
        self.__sensitivity_plot((3,2,3), 'sensitivity curve')

        self.__counts_plot(self.abs_all_counts,(3,2,4),'all positions',False,'upper left')
        self.__counts_plot(self.abs_ref_counts,(3,2,4),'matched positions',True,'upper left')

        top3_con,top3_val,top35_con,top35_val = self.__con_top()
        self.__plot_top_con(top3_con,top3_val,top35_con,top35_val,(3,2,5))

        self.__write_data()

        ##generate ref counts
        if self.args.ref_data:
            self.__read_ref_data()
            self.__plot_ref_data()


    def generate_all_plots(self):
        pdfout = PdfPages(self.args.output)
        self.plot = plt.figure(figsize=(10,10))
        self.label = 'Mutationseq (Single Sample) portraits for '+self.args.variant_file.split('/')[-1].replace('.vcf','')

        self.plot.suptitle(self.label,fontsize = 9)

        logging.info('Generating plots for: '+self.args.variant_file.split('/')[-1])

        self.__generate_plots()

        self.plot.subplots_adjust(left=0.05, bottom=0.03, right=0.97, top=0.93, wspace=0.2, hspace=0.2)

        #self.plot.tight_layout()
        pdfout.savefig(self.plot)
        pdfout.close()
        self.__delete_gzip()
        logging.info('Plotting finished successfully')

    def __con_top(self):
        if max(self.cum_all_counts.values()) >= 3000000:
            return self.__con_top_3mil()
        else:
            self.con_35_ref = False
            return self.__con_top_perc()

    def __con_top_perc(self):
        try:
            concordance_top_90 = float(self.cum_ref_counts['0.90'])/self.cum_all_counts['0.90']
            concordance_top_85 = float(self.cum_ref_counts['0.85'])/self.cum_all_counts['0.85']
        except:
            concordance_top_90 = 0
            concordance_top_85 = 0
            logging.error('Couldn\'t retreive the concordance info for .85/.90')
        return concordance_top_85,'0.85%',concordance_top_90,'0.90%'

    def __con_top_3mil(self):
        pr = None
        #get the values closest to top 3 mill
        minval_3 = min([abs(val[1]-3000000) for val in self.cum_all_counts.items()])
        val_3 = minval_3+3000000
        for k,val in self.cum_all_counts.iteritems():
            if val == val_3:
                pr = k

        if pr == None:
            val_3 = 3000000-minval_3
            for k,val in self.cum_all_counts.iteritems():
                if val == val_3:
                    pr = k

        concordance_3 = self.museq_percentages.get(pr)

        pr = None
        #get the values closest to top 3.5 mill
        minval_35 = min([abs(val[1]-3500000) for val in self.cum_all_counts.items()])
        val_35 = minval_35+3500000
        for k,val in self.cum_all_counts.iteritems():
            if val == val_35:
                pr = k

        if pr == None:
            val_35 = 3500000-minval_35
            for k,val in self.cum_all_counts.iteritems():
                if val == val_35:
                    pr = k

        concordance_35 = self.museq_percentages.get(pr)
        return concordance_3,val_3,concordance_35,val_35

    def __plot_top_con(self,top3_con,top3_val,top35_con,top35_val,subplot_axes,ref = False):
        logging.info('Generating concordance plot ')
        label_3 = 'top '+ str(top3_val)+' positions'
        label_35 = 'top '+ str(top35_val)+' positions'

        ax = self.plot.add_subplot(subplot_axes[0],subplot_axes[1],subplot_axes[2])

        if ref:
            ax.bar(1,top3_con,color=self.colors[1], label= 'reference '+label_3 )
            ax.bar(3,top35_con,color=self.colors[3], label= 'reference '+label_35 )
        else:
            ax.bar(0,top3_con,color=self.colors[0], label= label_3 )
            ax.bar(2,top35_con,color=self.colors[2], label= label_35 )

        if self.con_35_ref:
            ax.set_title(self.ref_labels+' Concordance for top 3/3.5 million positions', fontsize=8)
        else:
            ax.set_title(self.ref_labels+' Concordance for positions with probability 0.85/0.90 or above', fontsize=8)

        ax.set_ylabel("Concordance", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)

        ax.xaxis.set_tick_params(labelsize=0)
        ax.legend(loc = 'lower right',
                  fontsize = 6)


    def __concordance_plot(self,museq_percentages,subplot_axes,label):
        logging.info('Generating concordance plot ')
        title_con = self.ref_labels +' Concordance plot '
        keys = sorted(museq_percentages.keys())
        ax = self.plot.add_subplot(subplot_axes[0],subplot_axes[1],subplot_axes[2])
        ax.set_title(title_con, fontsize=8)
        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("percentage", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.plot(keys, [museq_percentages[k] for k in keys], linewidth=1.5,label = label)
        plt.ylim(0,100)
        ax.legend(loc = 'lower right',fontsize = 6,handlelength = 5)


    def __counts_plot(self,counts,subplot_axes,label,flag,loc='upper right'):
        logging.info('generate the counts_plot')
        if subplot_axes == (3,2,4):
            title_count = 'Density of all (solid) and matched (Dashed) positions in '+ self.ref_labels
        else:
            title_count = 'Number of snvs for all positions (solid) and \n positions found (Dashed) in '+ self.ref_labels
        keys = sorted(counts.keys())
        ax = self.plot.add_subplot(subplot_axes[0],subplot_axes[1],subplot_axes[2])
        ax.set_title(title_count, fontsize=8)
        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)
        if flag:
            ax.plot(keys, [counts[k] for k in keys],linestyle='--',label = label)
        else:
            ax.plot(keys, [counts[k] for k in keys], label = label)
        ax.legend(loc = loc,fontsize = 4,handlelength = 5)

    def __sensitivity_plot(self,subplot_axes,label,sensitivity= None):
        #matched_pos/all_pos in
        counts_all = self.len_ref
        logging.info('Calculating the sensitivity')
        title_sens = 'Sensitivity - percentage of '+self.ref_labels+' positions found'

        #calculate sensitivity if not provided
        if sensitivity == None:
            keys_ref = sorted(set(self.cum_ref_counts.keys()) )
            sensitivity = []
            for k in keys_ref:
                sensitivity.append(float(self.cum_ref_counts.get(k))/counts_all)
                self.sensitivity[k] = float(self.cum_ref_counts.get(k))/counts_all
        else:
            keys_ref = sorted(sensitivity.keys())
            sensitivity = [sensitivity[k] for k in keys_ref]


        ax = self.plot.add_subplot(subplot_axes[0],subplot_axes[1],subplot_axes[2])
        ax.set_title(title_sens, fontsize=8)
        ax.set_xlabel("probability", fontsize=8)
        ax.set_ylabel("Sensitivity(%)", fontsize=8)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.xaxis.set_tick_params(labelsize=6)
        ax.plot(list(keys_ref), sensitivity,linewidth=1.5,label = label)
        ax.legend(loc = 'lower left',fontsize = 6,handlelength = 5)

    # write all the data to a file
    def __write_data(self):
        self.outfile.write('\ncounts for:'+self.label+'\n')
        keys = sorted(self.cum_all_counts.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.cum_all_counts.get(k))+'\n')

        self.outfile.write('\nmatched counts for:'+self.label+'\n')
        keys = sorted(self.cum_ref_counts.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.cum_ref_counts.get(k))+'\n')

        self.outfile.write('\ndensity for:'+self.label+'\n')
        keys = sorted(self.abs_all_counts.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.abs_all_counts.get(k))+'\n')

        self.outfile.write('\nmatched density for:'+self.label+'\n')
        keys = sorted(self.abs_ref_counts.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.abs_ref_counts.get(k))+'\n')


        self.outfile.write('\nconcordance percentages for:'+self.label+'\n')
        keys = sorted(self.museq_percentages.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.museq_percentages.get(k))+'\n')

        self.outfile.write('\nsensitivity percentages for:'+self.label+'\n')
        keys = sorted(self.sensitivity.keys())
        for k in keys:
            self.outfile.write(str(k)+'\t'+str(self.sensitivity.get(k))+'\n')

    def __delete_gzip(self):
        try:
            for val in self.delfiles:
                os.remove(val)
        except:
            logging.error('Couldn\'t delete file: '+str(val))
