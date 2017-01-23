#!/usr/bin/env python

import os,sys
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
parser.add_argument("-s","--statsreg",type=str,
                    help="image region to calculate stats on")
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

imstat = BokImStat(extensions=extns,quickprocess=args.oscan,
                   stats_region=args.statsreg,checkbad=args.checkbad)
imstat.process_files(args.inputFiles)

#np.set_printoptions(precision=2,suppress=True)
np.set_printoptions(threshold=np.inf)

printarr = lambda a,fmt: ' '.join(map(lambda x: fmt%x,a))

for i,fn in enumerate(args.inputFiles):
	print os.path.basename(fn) 
	meanv = imstat.meanVals[i].reshape(4,4)
	rmsv = imstat.rmsVals[i].reshape(4,4)
	for j in range(4):
		print '  ',printarr(meanv[j],'%8.2f'),
		print '  ',printarr(rmsv[j],'%8.2f')
	if args.checkbad:
		nbad = np.array([badv.size for badv in imstat.badVals[i]]).reshape(4,4)
		for j in range(4):
			print '  ',printarr(nbad[j],'%8d')
		if args.showbad:
			for ext,badv in enumerate(imstat.badVals[i],start=1):
				print 'HDU[%d]:' % ext
				print badv.transpose()


