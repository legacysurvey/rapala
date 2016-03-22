#!/usr/bin/env python

import os,sys
import numpy as np
from astropy.stats import sigma_clip
import fitsio
from bokpipe import bokutil

class BokImStat(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokImStat,self).__init__(**kwargs)
		self.clipArgs = kwargs.get('clip_args',{})
	def _preprocess(self,fits,f):
		fn = os.path.basename(f)
		print '\n%s  ' % fn,
	def process_hdu(self,extName,data,hdr):
		pix = sigma_clip(data,**self.clipArgs)
		print '%8.2f ' % pix.mean(),
		return data,hdr

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inputFiles",type=str,nargs='+',
                    help="input FITS images")
parser.add_argument("-e","--ext",type=str,default="all",
                    help="select FITS extension(s) [default=all]")
args = parser.parse_args()

extns = args.ext.split(',')

imstat = BokImStat(extensions=extns)
imstat.process_files(args.inputFiles)


