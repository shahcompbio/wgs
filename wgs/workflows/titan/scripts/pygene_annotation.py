# Author:  Celia Siu
# Created: 13/11/2013

import pygene_annotation_ui
import sys
import os
import pygenes

args = pygene_annotation_ui.args

version = '1.0.2'

class pygeneAnnotation(object):
    def __init__(self,args):
        self.args = args
        ##checks: infile    
        if (not os.path.isfile(self.args.infile)):
            sys.stderr.write("error: '%s' is not a valid file\n" % args.infile)
            sys.exit()
        #checks: outfile
        if (self.args.outfile == None):
            self.outfile = os.path.abspath(self.args.infile) + ".pygenes"
        else:
            self.outfile = os.path.abspath(self.args.outfile)
    
        if (not os.path.exists(os.path.dirname(self.outfile))):
            os.makedirs(os.path.dirname(self.outfile))
            
        self.outfile = open(self.outfile, 'w')

    def save_bin(self):
        if self.args.gtf_bin:
            sys.stderr.write("Error: Can't save since the binary provided at input")
            sys.exit()
        gene_models = pygenes.GeneModels()
        gene_models.load_ensembl_gtf(args.gtf)
        gene_models.save_binary(self.args.outfile + '.gene_models.binary')
    
    def write_output(self):
        gene_models = pygenes.GeneModels()

        if (self.args.gtf_bin):
            gene_models.load_binary(self.args.gtf_bin)
        else:
            gene_models.load_ensembl_gtf(self.args.gtf)

        with open(args.infile, 'r') as titan_output:
            while True:
                line = titan_output.readline()
                if line[0] == '#':
                    continue
                header = line.rstrip()
                self.outfile.write("%s\tPygenes(gene_id,gene_name;)\n" % header)
                break

            ##do something to.. data lines
            for row in titan_output:
                row = row.rstrip()
                col = row.split('\t')

                if self.args.demix:
                    chrom = col[2]
                else:
                    chrom = col[1]
                try:
                    if self.args.demix:
                        start = int(col[3])
                        end = int(col[4])
                    else:
                        start = int(col[2])
                        end = int(col[3])
                except:
                    self.outfile.write("%s\t%s\n" % (row, "[ERROR - position not of type int()]"))
                    continue

                if (args.is_contained):
                    gene_ids = gene_models.find_contained_genes(chrom, start, end)
                else:
                    gene_ids = gene_models.find_overlapping_genes(chrom, start, end)

                pygenes_addition = ""
                for gene_id in gene_ids:
                    gene_name = gene_models.get_gene(gene_id).name
            
                    pygenes_addition += "%s,%s;" % (gene_id, gene_name)
                self.outfile.write("%s\t%s\n" % (row, pygenes_addition))
            self.outfile.close()
            
            
def _main():
    pyann = pygeneAnnotation(args)
    if args.gtf_save:
        pyann.save_bin()
    else:
        pyann.write_output()
        
if __name__ == '__main__':
    _main()
