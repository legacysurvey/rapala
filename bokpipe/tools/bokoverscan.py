#!/usr/bin/env/python

import sys
import argparse

from bokpipe.bokio import FileNameMap
from bokpipe.bokoscan import BokOverscanSubtract,BokOverscanSubtractWithSatFix

parser = argparse.ArgumentParser()
parser.add_argument("images",type=str,nargs='+',
                    help="input FITS images")
parser.add_argument("-f","--apply_filter",
                    help="([median]|none")
parser.add_argument("-m","--method",type=str,
                    help="([mean]|mean_value|median_value|cubic_spline)")
parser.add_argument("-r","--reject",type=str,
                    help="([sigma_clip]|minmax)")
parser.add_argument("--filter_window",type=int,
                    help="size of window to use with apply_filter")
parser.add_argument("--spline_nknots",type=int,
                    help="number of knots in spline fitting")
parser.add_argument("--spline_niter",type=int,
                    help="number of iterations of spline fitting")
parser.add_argument("--write_overscan_image",action="store_true",
                    help="make overscan image")
parser.add_argument("--oscan_cols_file",type=str,default="oscan_cols",
                    help="")
parser.add_argument("--oscan_rows_file",type=str,default="oscan_rows",
                    help="")
parser.add_argument("--row_apply_filter",
                    help="([median]|none")
parser.add_argument("--row_method",type=str,
                    help="([mean]|mean_value|median_value|cubic_spline)")
parser.add_argument("--row_reject",type=str,
                    help="([sigma_clip]|minmax)")
parser.add_argument("--row_filter_window",type=int,
                    help="size of window to use with apply_filter")
parser.add_argument("--row_spline_nknots",type=int,
                    help="number of knots in spline fitting")
parser.add_argument("--row_spline_niter",type=int,
                    help="number of iterations of spline fitting")
parser.add_argument('-p','--processes',type=int,default=1,
                    help='number of processes to use [default=single]')
args = parser.parse_args()

opts = vars(args)
kwargs = { k : opts[k] for k in opts if opts[k] != None }

omap = FileNameMap(newDir='.',newSuffix='_o')

oscan = BokOverscanSubtract(output_map=omap,**kwargs)

oscan.process_files(args.images)


