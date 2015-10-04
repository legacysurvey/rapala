#!/usr/bin/env python

import sys
import numpy as np
from astropy.stats import sigma_clip
import fitsio
import bokutil

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

imstat = BokImStat()

imstat.process_files(sys.argv[1:])


