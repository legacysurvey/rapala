#!/usr/bin/env/python

import sys
import argparse

from bokpipe.bokastrom import scamp_solve

parser = argparse.ArgumentParser()
parser.add_argument("image",type=str,
                    help="input FITS image")
parser.add_argument("catalog",type=str,
                    help="input FITS catalog")
parser.add_argument("-f","--filter",type=str,default='g',
                    help="reference band")
parser.add_argument("-o","--output",type=str,default=None,
                    help="output WCS header")
parser.add_argument("-p","--plots",action="store_true",
                    help="write check plots")
parser.add_argument("-r","--reference",type=str,default=None,
                    help="reference catalog")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
parser.add_argument("-w","--write",action="store_true",
                    help="write WCS to image header")
args = parser.parse_args()

if args.output is None:
	wcsFile = args.catalog.replace('.fits','.ahead')
else:
	wcsFile = args.output

scamp_solve(args.image,args.catalog,refStarCatFile=args.reference,
            filt=args.filter,savewcs=args.write,overwrite=True,
            check_plots=args.plots,verbose=args.verbose)

