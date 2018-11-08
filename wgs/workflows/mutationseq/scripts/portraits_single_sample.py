'''
Created on Mar 21, 2014

@author: dgrewal
@last modified: 3 Feb 2015 by jrosner
'''
import logging

import portraits_single_sample_utils
import portraits_single_sample_ui as ui

version = '1.0.1'

args = ui.args

if args.verbose:
    level = logging.DEBUG

else:
    level = logging.WARNING

logging.basicConfig(filename=args.log_file,
                    format='%(asctime)s %(message)s',
                    level=level)

logging.warning("<<< mutationSeq single sample concordance plotting started >>>")
logging.info(args)

concordance_plot = portraits_single_sample_utils.concordance_plots(args)
concordance_plot.generate_all_plots()
