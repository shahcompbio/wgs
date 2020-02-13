#!/share/lustre/csiu/apps/python2.7/bin/python2.7
# Author:  Celia Siu
# Created: 13/11/2013

import argparse

usage = """ %s [options] -i INFILE
        Note
        ----
        (I)
        - default HG19 gene sets were ftp transferred from:
          'ftp://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz'
        - ftp transfered file is found at:
          '/share/lustre/backup/scripts/pygenes/Homo_sapiens.GRCh37.73.gtf'
        - pygenes friendly binary format of this file is found at:
          '/share/lustre/backup/scripts/pygenes/Homo_sapiens.GRCh37.73.gtf.gene_models.binary'
        
        (II)
        - To save your own pygenes friendly gene sets binary format file, use the '-s' flag
          %s -i INFILE -r GENE_SETS.gtf -s
          You can use any INFILE that exists, this file is NOT read
          The binary output will be 'GENE_SETS.gtf.gene_models.binary'

        (III)
        - Adding the '-c' flag will return only genes that are contained 100 percent within segment;
          else results will include genes overlapping non-segment regions at the start and ends.
        """

parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--infile', 
                    dest='infile',
                    required=True,
                    help='path to input file ... the output of Titan')

parser.add_argument('-o', '--outfile', 
                    dest='outfile',
                    help="output file name; Default = INFILE.pygenes")

parser.add_argument('-s', '--save_gtf_as_bin', 
                    dest='gtf_save',
                    action='store_true',
                    help="save '-r' file to pygenes friendly binary format")

parser.add_argument('-c', '--is_contained', 
                    dest='is_contained',
                    action='store_true',
                    help='add flag to return output that is contained 100 percent within region \n\
                         (else will include genes overlapping segment ends)')

parser.add_argument( '--demix',
                    dest='demix',
                    action='store_true',
                    default = False,
                    help='add flag to set input file format to demix')

gtf_infile_group = parser.add_mutually_exclusive_group(required = True)
gtf_infile_group.add_argument('-b', '--gene_sets_gtf_bin',
                              dest='gtf_bin',
                              help="Gene sets in GTF binary format")   
                 
gtf_infile_group.add_argument('-r', '--gene_sets_gtf',
                              dest='gtf',
                              help="Gene sets in GTF format")

##get at the arguments
args = parser.parse_args()
