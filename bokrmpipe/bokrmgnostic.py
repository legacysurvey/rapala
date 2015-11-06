#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import fitsio

from bokpipe import *

def plot_gain_vals(diagfile):
	g = np.load(diagfile)#,gains=gainCorV,skys=skyV,gainCor=gainCor)
	plt.figure(figsize=(9,6))
	plt.subplots_adjust(0.07,0.04,0.97,0.97,0.25,0.05)
	for amp in range(16):
		ax = plt.subplot(4,4,amp+1)
		plt.plot(g['gains'][:,0,amp],c='b')
		plt.axhline(g['gainCor'][0,amp],c='purple',ls='--')
		plt.plot(g['gains'][:,1,amp],c='r')
		plt.axhline(g['gainCor'][1,amp],c='orange',ls='--')
		ax.xaxis.set_visible(False)
		plt.ylim(0.91,1.09)
		ax.text(0.05,0.05,'IM%d'%bokproc.ampOrder[amp],
		        size=8,transform=ax.transAxes)
		ax.text(0.25,0.05,'%.3f'%g['gainCor'][0,amp],color='blue',
		        size=8,transform=ax.transAxes)
		ax.text(0.50,0.05,'%.3f'%g['gainCor'][1,amp],color='red',
		        size=8,transform=ax.transAxes)

from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage.filters import maximum_filter

def make_circular_filter(radius):
	y,x = np.indices((radius*2,radius*2),dtype=np.float32)
	r = np.sqrt((x+0.5-radius)**2 + (y+0.5-radius)**2)
	circfilt = r < radius
	return circfilt

def make_gradient_ims(files,bpmask,nbin=16,nims=None,
                      sfx='',extensions=None,dorms=False):
	if extensions is None:
		extensions = ['IM9','IM10']
	binnedIms = { extn:[] for extn in extensions }
	snum = 1
	if nims is None:
		nims = len(files)
	for f in files:
		print 'processing ',f
		mef = bokutil.BokMefImage(f,read_only=True,mask_file=bpmask)
		outf = 'stack%s_%d.fits' % (sfx,snum)
		if os.path.exists(outf):
			os.unlink(outf)
		outFits = fitsio.FITS(outf,'rw')
		outFits.write(None)
		if dorms:
			outrmsf = 'stackrms%s_%d.fits' % (sfx,snum)
			if os.path.exists(outrmsf):
				os.unlink(outrmsf)
			outrmsFits = fitsio.FITS(outrmsf,'rw')
			outrmsFits.write(None)
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
			if len(binnedIms[extn])==nims:
				cube = np.ma.dstack(binnedIms[extn])
				cube = sigma_clip(cube,iters=3,sig=2.0,axis=-1)
				stack = cube.mean(axis=-1).filled(0)
				stackrms = cube.std(axis=-1).filled(0)
				print 'writing stack ',snum
				outFits.write(stack,extname=extn,header={'NBIN':nbin})
				if dorms:
					outrmsFits.write(stackrms,extname=extn)
				if extn==extensions[-1]:
					snum += 1
				binnedIms[extn] = []
		outFits.close()
		if dorms:
			outrmsFits.close()

def make_correction_im(gradientFile,outputFile,fitFile=None):
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
		ii = gim != 0
		spfit = LSQBivariateSpline(xi[ii]*nbin,yi[ii]*nbin,gim[ii],
		                           [0,nX],[0,nY],kx=3,ky=3)
		im = spfit(np.arange(nX),np.arange(nY)).T
		corrFits.write(im,extname=extn)
		if fitFile is not None:
			fitFits.write(spfit(np.arange(nx)*nbin,np.arange(ny)*nbin).T,
			              extname=extn)
	gradientFits.close()
	corrFits.close()
	if fitFile is not None:
		fitFits.close()

def make_darksky_correction():
	# XXX need to actually pull correct files from log
	from glob import glob
	pdir='tmprm_bias/'
	f0427 = sorted(glob(pdir+'ut20140427/bokrm.20140427.????.fits'))
	f0427 = f0427[10:]
	bpmaskf = '/d1/OldDataDisk1/dev/rapala/bokrmpipe/tmprm/cals/badpix_master.fits'
	make_gradient_ims(f0427,bpmaskf,sfx='_p0427')
	gradientFile = pdir+'diagnostics/stack_p0427_1.fits'
	outputFile = pdir+'cals/biasramp.fits'
	make_correction_im(gradientFile,outputFile)#,fitFile=None)

def test_gradient_ims(pdir='tmprm_bias/'):
	from glob import glob
	f0123 = sorted(glob(pdir+'ut20140123/bokrm.20140123.????.fits'))
	f0123 = f0123[27:57]
	make_gradient_ims(f0123,sfx='_p0123')
	f0413 = sorted(glob(pdir+'ut20140413/bokrm.20140413.????.fits'))
	f0413 = f0413[30:60]
	make_gradient_ims(f0413,sfx='_p0413')
	f0415 = sorted(glob(pdir+'ut20140415/bokrm.20140415.????.fits'))
	f0415 = f0415[12:43]
	make_gradient_ims(f0415,sfx='_p0415')
	f0427 = sorted(glob(pdir+'ut20140427/bokrm.20140427.????.fits'))
	f0427 = f0427[10:]
	make_gradient_ims(f0427,sfx='_p0427')

if __name__=='__main__':
	make_darksky_correction()

