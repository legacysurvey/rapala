#!/usr/bin/env python

import sys
import numpy as np
import fitsio
import bokutil

class BokDiff(bokutil.BokProcess):
	def __init__(self,file1,file2,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokDiff,self).__init__(**kwargs)
		self.fits2 = fitsio.FITS(file2)
		self.isSame = True
		self.file1,self.file2 = file1,file2
	def process_hdu(self,extName,data,hdr):
		self.isSame &= np.all(np.isclose(data,self.fits2[extName].read()))
		return data,hdr
	def _postprocess(self):
		if self.isSame:
			print '%s and %s are the same' % (self.file1,self.file2)
		else:
			print '%s and %s are different' % (self.file1,self.file2)

imdiff = BokDiff(*sys.argv[1:])

imdiff.process_files([sys.argv[1]])


