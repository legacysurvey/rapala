#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import fitsio

from bokpipe import *

from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage.filters import maximum_filter

def make_circular_filter(radius):
	y,x = np.indices((radius*2,radius*2),dtype=np.float32)
	r = np.sqrt((x+0.5-radius)**2 + (y+0.5-radius)**2)
	circfilt = r < radius
	return circfilt

def make_gradient_ims(files,outputFile,bpmask,nbin=16,
                      extensions=None,dorms=False):
	if extensions is None:
		extensions = ['IM9','IM10']
	binnedIms = { extn:[] for extn in extensions }
	if os.path.exists(outputFile):
		os.unlink(outputFile)
	outFits = fitsio.FITS(outputFile,'rw')
	outFits.write(None)
	if dorms:
		outrmsf = outputFile.replace('.fits','_rms.fits')
		if os.path.exists(outrmsf):
			os.unlink(outrmsf)
		outrmsFits = fitsio.FITS(outrmsf,'rw')
		outrmsFits.write(None)
	for f in files:
		print 'constructing ramp image for ',f
		mef = bokutil.BokMefImage(f,read_only=True,mask_file=bpmask)
		for extn in extensions:
			im = mef.get(extn)
			# do a global image clipping 
			im = sigma_clip(im,iters=3,sig=5.0)
			subim = im.data - im.mean()
			# construct masks around bright stars
			findmax = maximum_filter(subim,(150,150))
			y,x = np.indices(im.shape,dtype=np.float32)
			masked = np.zeros(im.shape,dtype=bool)
			for cut,rad in zip([20000.,10000.,5000.,2500],[50,30,15,10]):
				ii = np.where((subim==findmax) & (subim > cut) & ~masked)
				for _i,_j in zip(*ii):
					masked |= np.sqrt((x-_i)**2+(y-_j)**2) < rad
			im.mask |= masked
			# grow the mask a bunch
			im.mask |= binary_dilation(im.mask,iterations=10)
			# rebin the masked image
			bim = bokutil.rebin(im,nbin)
			# count number of masked pixels in each binned pixel
			nmasked = bim.mask.sum(axis=-1)
			# now do a local clipping within each bin
			bim = sigma_clip(bim,iters=3,sig=2.2,axis=-1)
			mbim = bim.mean(axis=-1)
			# mask if a majority of binned pixels are masked
			mbim.mask |= nmasked > nbin**2/2
			# subtract the median sky level
			mbim -= np.ma.median(mbim)
			binnedIms[extn].append(mbim)
	#
	for extn in extensions:
		cube = np.ma.dstack(binnedIms[extn])
		cube = sigma_clip(cube,iters=3,sig=2.0,axis=-1)
		stack = cube.mean(axis=-1).filled(0)
		stackrms = cube.std(axis=-1).filled(0)
		outFits.write(stack,extname=extn,header={'NBIN':nbin})
		if dorms:
			outrmsFits.write(stackrms,extname=extn)
	outFits.close()
	if dorms:
		outrmsFits.close()

def make_correction_im(gradientFile,outputFile,fitFile=None,nKnots=5):
	from scipy.interpolate import LSQBivariateSpline
	gradientFits = fitsio.FITS(gradientFile)
	if os.path.exists(outputFile):
		os.unlink(outputFile)
	corrFits = fitsio.FITS(outputFile,'rw')
	corrFits.write(None)
	if fitFile is not None:
		fitFits = fitsio.FITS(fitFile,'rw')
		fitFits.write(None)
	for hdu in gradientFits[1:]:
		extn = hdu.get_extname().upper()
		gim = hdu.read()
		nbin = hdu.read_header()['NBIN']
		ny,nx = gim.shape
		nY,nX = ny*nbin,nx*nbin
		yi,xi = np.indices(gim.shape)
		mask = gim == 0
		if extn=='IM9':
			# add extra mask around a group of hot pixels
			mask |= ( (xi*nbin > 190) & (xi*nbin < 290) &
			          (yi*nbin > 235) & (xi*nbin < 260) )
		ii = np.where(~mask)
		xknots = np.linspace(0,nX,nKnots)
		yknots = np.linspace(0,nY,nKnots)
		spfit = LSQBivariateSpline(xi[ii]*nbin,yi[ii]*nbin,gim[ii],
		                           xknots,yknots,kx=3,ky=3)
		im = spfit(np.arange(nX),np.arange(nY)).T
		corrFits.write(im,extname=extn)
		if fitFile is not None:
			fitFits.write(spfit(np.arange(nx)*nbin,np.arange(ny)*nbin).T,
			              extname=extn)
	gradientFits.close()
	corrFits.close()
	if fitFile is not None:
		fitFits.close()

def make_rampcorr_image(dataMap,**kwargs):
	gradientFile = os.path.join(dataMap._tmpDir,'tmpramp.fits')
	files = [dataMap('proc1',False)(f) for f in dataMap.getFiles()]
	if not os.path.exists(gradientFile):
		make_gradient_ims(files,gradientFile,dataMap('MasterBadPixMask'))
	outFn = os.path.join(dataMap.calDir,dataMap.masterRampCorrFn)
	make_correction_im(gradientFile,outFn)

