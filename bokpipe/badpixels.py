#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.morphology import binary_dilation,binary_closing
from astropy.stats import sigma_clip

from bokutil import rebin,BokMefImage

import fitsio


def build_mask_from_flat(flatFn,outFn,**kwargs):
	loCut,hiCut = kwargs.get('good_range',(0.92,1.08))
	loCut2,hiCut2 = kwargs.get('grow_range',(0.97,1.03))
	nbin = kwargs.get('nbin',32)
	flatName = kwargs.get('normed_flat_file')
	bpMask = BokMefImage(flatFn,output_file=outFn,**kwargs)
	if flatName is not None:
		if os.path.exists(flatName):
			os.unlink(flatName)
		normedFlat = fitsio.FITS(flatName,'rw')
		normedFlat.write(None,header=fitsio.read_header(flatFn,0))
	for extName,data,hdr in bpMask:
		ny,nx = data.shape
		x = np.arange(nbin/2,nx,nbin)
		y = np.arange(nbin/2,ny,nbin)
		im = np.ma.masked_array(data,mask=((data < 0.4) | (data>1.5)))
		binnedIm = rebin(im,nbin)
		binnedIm = sigma_clip(binnedIm,axis=-1,iters=3,sig=2.2,
		                      cenfunc=np.ma.mean).mean(axis=-1).filled(1)
		spfit = RectBivariateSpline(x,y,binnedIm.T,s=1)
		gradientIm = spfit(np.arange(nx),np.arange(ny)).T
		im /= gradientIm
		if flatName is not None:
			normedFlat.write(im,extname=extName,header=hdr)
		badpix = (im < loCut) | (im > hiCut)
		badpix |= binary_dilation(badpix,mask=((im<loCut2)|(im>hiCut2)),
		                          iterations=0)
		badpix |= binary_closing(badpix,iterations=0)
		bpMask.update(badpix.astype(np.int16),hdr)
	bpMask.close()
	if flatName is not None:
		normedFlat.close()

def make_sextractor_gain_map(flatFn,bpMaskFn,gainMapFn,**kwargs):
	flat = BokMefImage(flatFn,output_file=gainMapFn,**kwargs)
	mask = fitsio.FITS(bpMaskFn)
	for extName,data,hdr in flat:
		# invert the mask, makes bad pixels = 0
		data *= (1 - mask[extName].read())
		flat.update(data,hdr)
	flat.close()
