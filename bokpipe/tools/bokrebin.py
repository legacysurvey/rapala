#!/usr/bin/env python

import os,sys
import numpy as np
import fitsio
from bokpipe.bokutil import rebin,BokMefImage
from bokpipe.bokoscan import _convertfitsreg

fileName = sys.argv[1]
if len(sys.argv)==3:
	nbin = int(sys.argv[2])
	maskFile = None
elif len(sys.argv)==4:
	nbin,maskFile = int(sys.argv[2]),sys.argv[3]

fits = BokMefImage(fileName,mask_file=maskFile,read_only=True)

outFits = fitsio.FITS(fileName.replace('.fits','_bin%d.fits'%nbin),'rw',
                      clobber=True)

def convertreg(reg):
	a = np.array(_convertfitsreg(reg))
	reg = 1 + a//nbin
	return tuple(reg)

hdr = fits.get_header(0)
hdr['DETSIZE'] = '[%d:%d,%d:%d]' % convertreg(hdr['DETSIZE'])
outFits.write(None,header=hdr)

for extn,data,hdr in fits:
	im = rebin(data,nbin).astype(np.float32).mean(axis=-1)
	for c in ['CRPIX1','CRPIX2','LTV1','LTV2']:
		hdr[c] //= nbin
	for c in ['CD1_1','CD1_2','CD1_2','CD2_2','LTM1_1','LTM2_2']:
		hdr[c] *= nbin
	hdr['DETSEC'] = '[%d:%d,%d:%d]' % convertreg(hdr['DETSEC'])
	outFits.write(im,extname=extn,header=hdr)

outFits.close()


