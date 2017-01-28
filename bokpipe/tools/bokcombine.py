#!/usr/bin/env python

import os
import argparse

from bokpipe import *

def _read_file_list(files):
	rv = []
	# accept an IRAF-style list from a file if args starts with '@'
	for f in files:
		if f.startswith('@'):
			with open(f[1:]) as _f:
				rv.extend(_f.read().splitlines())
		else:
			rv.append(f)
	return rv

parser = argparse.ArgumentParser()
parser.add_argument("inputFiles",type=str,nargs='+',
                    help="input FITS images")
parser.add_argument("-o","--output",type=str,default="bokstack.fits",
                    help="output image")
parser.add_argument("-m","--method",type=str,default="mean",
                    help="stacking method (mean|median|mode)")
parser.add_argument("-s","--scale",action='store_true',
#                    default="normalize_median",
                    help="scale images before combining")
parser.add_argument("-r","--reject",action='store_true',
#                    default="sigma_clip",
                    help="apply rejection")
parser.add_argument("--maxmem",type=int,default=8,
                    help="maximum memory in GB")
parser.add_argument("-w","--weights",type=str,nargs='+',
                    help="weight images")
parser.add_argument("--masksfx",type=str,
                    help="input mask image sfx (e.g., '.msk')")
parser.add_argument("--maskdir",type=str,default='.',
                    help="directory where masks are located (default:.)")
parser.add_argument("--masktype",type=str,default='gtzero',
                    help="input mask boolean type ([gtzero]|nonzero|zerobad)")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
args = parser.parse_args()

stackPars = {}
if args.scale:
	stackPars['scale'] = 'normalize_mean'
if args.reject:
	stackPars['reject'] = 'sigma_clip'
stackPars['maxmem'] = args.maxmem
stackPars['verbose'] = args.verbose

if args.masksfx:
	stackPars['mask_map'] = lambda f: os.path.join(args.maskdir,
	                             f.replace('.fits',args.masksfx+'.fits'))
	stackPars['mask_type'] = args.masktype

if args.method == 'mean':
	stackFun = bokutil.ClippedMeanStack(**stackPars)
elif args.method == 'median':
	stackFun = bokutil.MedianStack(**stackPars)

inputFiles = _read_file_list(args.inputFiles)
print 'input: ',inputFiles
print 'output: ',args.output
print stackPars
if args.weights is not None:
	weightFiles = _read_file_list(args.weights)
	print 'weights: ',weightFiles
else:
	weightFiles = None

stackFun.stack(inputFiles,args.output,weights=weightFiles)

