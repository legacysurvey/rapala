#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.morphology import binary_dilation,binary_closing
from astropy.stats import sigma_clip

from bokutil import rebin,BokMefImage

def build_mask_from_flat(flatFn,outFn,**kwargs):
	loCut = kwargs.get('lo_cut',0.8)
	hiCut = kwargs.get('hi_cut',1.2)
	nbin = kwargs.get('nbin',32)
	bpMask = BokMefImage(flatFn,output_file=outFn,**kwargs)
	for data,hdr in bpMask:
		ny,nx = data.shape
		x = np.arange(nbin/2,nx,nbin)
		y = np.arange(nbin/2,ny,nbin)
		im = np.ma.masked_array(data,mask=((data < 0.4) & (data>1.5)))
		binnedIm = rebin(im,nbin)
		binnedIm = sigma_clip(binnedIm,axis=-1,iters=3,sig=2.2,
		                      cenfunc=np.ma.mean).mean(axis=-1).filled(0)
		spfit = RectBivariateSpline(x,y,binnedIm.T,s=1)
		gradientIm = spfit(np.arange(nx),np.arange(ny)).T
		im /= gradientIm
		#bpMask.update(im,hdr)
		badpix = (im < 0.8) | (im > 1.2)
		badpix |= binary_dilation(badpix,mask=((im<0.95)|(im>1.05)),
		                          iterations=0)
		badpix |= binary_closing(badpix)
		bpMask.update(badpix.astype(np.int16),hdr)
	bpMask.close()

