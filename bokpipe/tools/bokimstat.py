#!/usr/bin/env python

import os,sys
import numpy as np
from astropy.stats import sigma_clip
import fitsio

from bokpipe import bokutil
from bokpipe.bokoscan import overscan_subtract

class BokImStat(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokImStat,self).__init__(**kwargs)
		self.clipArgs = kwargs.get('clip_args',{})
		self.quickprocess = kwargs.get('quickprocess',False)
	def _preprocess(self,fits,f):
		fn = os.path.basename(f)
		print '\n%s  ' % fn,
	def process_hdu(self,extName,data,hdr):
		if self.quickprocess:
			pix = overscan_subtract(data,hdr,method='mean_value',
			                        reject='sigma_clip',clip_iters=1,
			                        apply_filter=None)
		else:
			pix = data
		pix = sigma_clip(pix,**self.clipArgs)
		print '%8.2f ' % pix.mean(),
		return data,hdr

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inputFiles",type=str,nargs='+',
                    help="input FITS images")
parser.add_argument("-e","--ext",type=str,default="all",
                    help="select FITS extension(s) [default=all]")
parser.add_argument("--oscan",action="store_true",
                    help="do a simple overscan subtraction")
args = parser.parse_args()

if args.ext=='all':
	extns = None
else:
	extns = args.ext.split(',')

imstat = BokImStat(extensions=extns,quickprocess=args.oscan)
imstat.process_files(args.inputFiles)


