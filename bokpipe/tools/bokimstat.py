#!/usr/bin/env python

import os,sys
import multiprocessing
import numpy as np
from astropy.stats import sigma_clip
import fitsio

from bokpipe.bokproc import BokImStat

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inputFiles",type=str,nargs='+',
                    help="input FITS images")
parser.add_argument("-e","--ext",type=str,default="all",
                    help="select FITS extension(s) [default=all]")
parser.add_argument("-f","--fields",type=str,default="mean,std",
                    help="which stats to calculate [default: mean,std]")
parser.add_argument("-s","--statsreg",type=str,
                    help="image region to calculate stats on")
parser.add_argument("--stride",type=int,
                    help="stride to use in stats region")
parser.add_argument("-p","--processes",type=int,default=1,
                    help="use multiple processes")
parser.add_argument("--oscan",action="store_true",
                    help="do a simple overscan subtraction")
parser.add_argument("--checkbad",action="store_true",
                    help="check for bad values in the image")
parser.add_argument("--showbad",action="store_true",
                    help="dump coordinates of bad values")
args = parser.parse_args()

if args.ext=='all':
	extns = None
else:
	extns = args.ext.split(',')

if args.processes > 1:
	pool = multiprocessing.Pool(args.processes)
	procmap = pool.map
else:
	procmap = map

fields = args.fields.split(',')

imstat = BokImStat(extensions=extns,fields=fields,
                   stats_region=args.statsreg,stats_stride=args.stride,
                   quickprocess=args.oscan,checkbad=args.checkbad,
                   processes=args.processes,procmap=procmap)
imstat.process_files(args.inputFiles)

#np.set_printoptions(precision=2,suppress=True)
np.set_printoptions(threshold=np.inf)

printarr = lambda a,fmt: ' '.join(map(lambda x: fmt%x,a))

if extns is None:
	nrow = imstat.data[fields[0]].shape[1] // 4
	ncol = 4
else:
	nrow = len(extns)
	ncol = 1

for i,fn in enumerate(args.inputFiles):
	print os.path.basename(fn),
	for k in fields:
		v = imstat.data[k][i].reshape(nrow,ncol)
		for j in range(nrow):
			print '  ',printarr(v[j],'%8.2f'),
	print
	if args.checkbad:
		nbad = np.array([badv.size for badv in imstat.badVals[i]]).reshape(nrow,ncol)
		for j in range(nrow):
			print '  ',printarr(nbad[j],'%8d')
		if args.showbad:
			for ext,badv in enumerate(imstat.badVals[i],start=1):
				print 'HDU[%d]:' % ext
				print badv.transpose()


