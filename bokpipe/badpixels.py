#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.morphology import binary_dilation,binary_closing
from astropy.stats import sigma_clip
import fitsio

from bokutil import rebin,BokMefImage,BokProcess

class BadPixelMaskFromFlats(BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','BPMSK')
		super(BadPixelMaskFromFlats,self).__init__(**kwargs)
		self.loCut,self.hiCut = kwargs.get('good_range',(0.90,1.10))
		self.loCut2,self.hiCut2 = kwargs.get('grow_range',(0.95,1.05))
		self.nbin = kwargs.get('nbin',32)
		self.flatName = kwargs.get('normed_flat_file')
		self.flatFitName = kwargs.get('normed_flat_fit_file')
		self.binnedFlatName = kwargs.get('binned_flat_file')
		self.noConvert = True  # casting to unsigned int below
		self.normedFlat = None
		self.normedFlatFit = None
		self.binnedFlat = None
	def _preprocess(self,fits,f):
		if self.flatName is not None:
			if os.path.exists(self.flatName):
				os.unlink(self.flatName)
			self.normedFlat = fitsio.FITS(self.flatName,'rw')
			self.normedFlat.write(None,header=fits.get_header(0))
		if self.flatFitName is not None:
			if os.path.exists(self.flatFitName):
				os.unlink(self.flatFitName)
			self.normedFlatFit = fitsio.FITS(self.flatFitName,'rw')
			self.normedFlatFit.write(None,header=fits.get_header(0))
		if self.binnedFlatName is not None:
			if os.path.exists(self.binnedFlatName):
				os.unlink(self.binnedFlatName)
			self.binnedFlat = fitsio.FITS(self.binnedFlatName,'rw')
			self.binnedFlat.write(None,header=fits.get_header(0))
	def _postprocess(self):
		if self.normedFlat is not None:
			self.normedFlat.close()
			self.normedFlat = None
		if self.normedFlatFit is not None:
			self.normedFlatFit.close()
			self.normedFlatFit = None
		if self.binnedFlat is not None:
			self.binnedFlat.close()
			self.binnedFlat = None
	def process_hdu(self,extName,data,hdr):
		ny,nx = data.shape
		im = np.ma.masked_array(data,mask=((data<0.5) | (data>1.5)))
		# mask the edges
		margin = 15
		im.mask[:margin] = True
		im.mask[-margin:] = True
		im.mask[:,:margin] = True
		im.mask[:,-margin:] = True
		if extName == 'IM9':
			# this one has a long bad strip
			im.mask[:,:50] = True
		binnedIm = rebin(im,self.nbin)
		binnedIm = sigma_clip(binnedIm,axis=-1,iters=3,sig=2.2,
		                      cenfunc=np.ma.mean).mean(axis=-1).filled(1)
		x = np.arange(self.nbin/2,nx,self.nbin)
		y = np.arange(self.nbin/2,ny,self.nbin)
		spfit = RectBivariateSpline(x,y,binnedIm.T,s=1)
		gradientIm = spfit(np.arange(nx),np.arange(ny)).T
		im /= gradientIm
		if self.normedFlat is not None:
			self.normedFlat.write(im.astype(np.float32),
			                      extname=extName,header=hdr)
		if self.normedFlatFit is not None:
			self.normedFlatFit.write(gradientIm.astype(np.float32),
			                         extname=extName,header=hdr)
		if self.binnedFlat is not None:
			self.binnedFlat.write(binnedIm.astype(np.float32),
			                      extname=extName,header=hdr)
		badpix = (im < self.loCut) | (im > self.hiCut)
		badpix |= binary_dilation(badpix,
		                          mask=((im<self.loCut2)|(im>self.hiCut2)),
		                          iterations=0)
		badpix |= binary_closing(badpix,iterations=0)
		return badpix.astype(np.uint8),hdr

def build_mask_from_flat(flatFn,outFn,**kwargs):
	bpmGen = BadPixelMaskFromFlats(output_map=lambda f: outFn,
	                               **kwargs)
	bpmGen.process_files([flatFn,])



def make_sextractor_gain_map(flatFn,bpMaskFn,gainMapFn,**kwargs):
	flat = BokMefImage(flatFn,output_file=gainMapFn,**kwargs)
	mask = fitsio.FITS(bpMaskFn)
	for extName,data,hdr in flat:
		# invert the mask, makes bad pixels = 0
		data *= (1 - mask[extName].read())
		flat.update(data,hdr)
	flat.close()
