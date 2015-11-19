#!/usr/bin/env python

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
parser.add_argument("-s","--scale",type=str,default="normalize_median",
                    help="scale images before combining")
parser.add_argument("-r","--reject",type=str,default="sigma_clip",
                    help="apply rejection")
parser.add_argument("-n","--nsplit",type=int,default=8,
                    help="number of chunks to split images into")
parser.add_argument("-w","--weights",type=str,nargs='+',
                    help="weight images")
args = parser.parse_args()

stackPars = {}
stackPars['scale'] = args.scale
stackPars['reject'] = args.reject
stackPars['nsplit'] = args.nsplit

if args.method == 'mean':
	stackFun = bokutil.ClippedMeanStack(**stackPars)
elif args.method == 'median':
	stackFun = bokutil.MedianStack(**stackPars)

inputFiles = _read_file_list(args.inputFiles)
print 'input: ',inputFiles
print 'output: ',args.output
if args.weights is not None:
	weightFiles = _read_file_list(args.weights)
	print 'weights: ',weightFiles

stackFun.stack(inputFiles,args.output,weights=weightFiles)

