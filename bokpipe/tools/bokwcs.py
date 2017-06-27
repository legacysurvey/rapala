#!/usr/bin/env/python

import sys
import argparse

from bokpipe.bokastrom import scamp_solve

parser = argparse.ArgumentParser()
parser.add_argument("image",type=str,
                    help="input FITS image")
parser.add_argument("catalog",type=str,
                    help="input FITS catalog")
parser.add_argument("-a","--args",type=str,
                    help="arguments to pass to scamp config")
parser.add_argument("-f","--filter",type=str,default='g',
                    help="reference band")
parser.add_argument("-p","--plots",action="store_true",
                    help="write check plots")
parser.add_argument("-r","--reference",type=str,default=None,
                    help="reference catalog")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
parser.add_argument("-w","--write",action="store_true",
                    help="write WCS to image header")
parser.add_argument("--single",action="store_true",
                    help="single pass")
args = parser.parse_args()

kwargs = {}
if args.args is not None:
	arglist = args.args.split()
	for a in arglist:
		k,v = a.split('=')
		kwargs[k] = v

scamp_solve(args.image,args.catalog,refStarCatFile=args.reference,
            filt=args.filter,savewcs=args.write,clobber=True,
            check_plots=args.plots,twopass=not args.single,
            verbose=args.verbose,**kwargs)

