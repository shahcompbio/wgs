'''
Created on Mar 20, 2014

@author: dgrewal
@last modified: 3 Feb 2015 by jrosner
'''

import argparse

parser = argparse.ArgumentParser(prog='generate feature boxplots and concordance plots',
                                 description='''generates boxplots for feature
                                                distributions and concordance plots''')

parser.add_argument('-v','--variant_file',
                    required = True,
                    help = 'Path to the variant file (/path/to/vcf)')

parser.add_argument('-b','--variant_label',
                    required = True,
                    help = 'The label for the variant. Used to label the axis on the plots.')

parser.add_argument("-l", "--log_file",
                    default="portraits_single_sample_run.log",
                    help='''specify name or path of the log file''')

parser.add_argument('--output',
                    required = True,
                    help = '''path to the folder with prefix for the
                            outputs (/path/to/output_folder/name)''')

parser.add_argument('--data',
                    required = True,
                    help = '''path to the folder with prefix for the
                            outputs (/path/to/output_folder/name)''')

parser.add_argument('--dbsnp',
                    help='path to dbsnp reference file')

parser.add_argument('--thousand_gen',
                    help='path to thousand genomes reference file')

parser.add_argument('--tabix_path',
                    default = None,
                    help = '''path to the tabix installation folder,
                              if the inputs are vcf then they will be indexed by tabix ''')

parser.add_argument('--threshold',
                    default = 0.5,
                    type = float,
                    help = ''' threshold''')

parser.add_argument( "--verbose",
                    action="store_true",
                    default=False,
                    help='''verbose''')

parser.add_argument("-d","--ref_data",
                    default = None,
                    help = '''Adds reference information to plots''')

args = parser.parse_args()
