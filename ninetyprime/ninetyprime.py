#!/usr/bin/env python

import re
import numpy as np
from astropy.stats import sigma_clip

nX,nY = 4032,4096

# the order of the amplifiers in the FITS extensions, i.e., HDU1=amp#4
ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]


# extremely basic CCD image processing routines for quick look purposes

def _convertfitsreg(regstr):
	regpattern = r'\[(\d+):(\d+),(\d+):(\d+)\]'
	rv =  [ int(d) for d in  re.match(regpattern,regstr).groups() ]
	# FITS region indices are 1-indexed
	rv[0] -= 1
	rv[2] -= 1
	return rv

def colbias(imhdu,method='median',xmargin1=5,xmargin2=2,ymargin=10,**kwargs):
	data = imhdu.read()
	hdr = imhdu.read_header()
	x1,x2,y1,y2 = _convertfitsreg(hdr['DATASEC'])
	im = data[y1:y2,x1:x2].astype(np.float32)
	x1,x2,y1,y2 = _convertfitsreg(hdr['BIASSEC'])
	bias = data[y1:y2,x1:x2].astype(np.float32)
	bias = sigma_clip(bias[ymargin:-ymargin,xmargin1:-xmargin2])
	# different versions of astropy?
	try:
		bias = np.ma.mean(bias).filled()[0]
	except:
		bias = np.ma.mean(bias)
	im -= bias
	return im,bias

def improcess(imhdu,extn=None,biasim=None,pixflatim=None,superskyim=None,
              **kwargs):
	im,bias = colbias(imhdu,**kwargs)
	if biasim is not None:
		im -= biasim[extn]
	if pixflatim is not None:
		im /= pixflatim[extn].data
	if superskyim is not None:
		im /= superskyim[extn]
	return im,bias

