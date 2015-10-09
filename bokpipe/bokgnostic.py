#!/usr/bin/env python

import numpy as np
import fitsio

import matplotlib.pyplot as plt

import bokutil
import bokproc

def check_gain_bal(fileName,badPixMaskFile=None,
                   showMean=True,showMode=True,**kwargs):
	maskFits = None
	if badPixMaskFile is not None:
		maskFits = fitsio.FITS(badPixMaskFile)
	gainBalance = bokproc.BokCalcGainBalanceFactors(save_arrays=True,
                                                 mask_map=lambda f: maskFits,
	                                                **kwargs)
	gainBalance.process_files([fileName])
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(16):
		ax = plt.subplot(4,4,i+1)
		extn = 'IM%d'%(i+1)
		pix = gainBalance.arrays[i].flatten()
		modalVal = 3*np.ma.median(pix)-2*pix.mean()
		plt.hist(pix,100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(gainBalance.skyVals[i],color='r',lw=1.2)
		if showMean:
			plt.axvline(pix.mean(),color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,
		         np.ma.median(pix),pix.mean(),
		         pix.std(),pix.min(),pix.max())

def check_sky_level(fileName,maskFile=None,statreg='ccd_central_quadrant',
                    showMean=True,showMode=True,**kwargs):
	fits = fitsio.FITS(fileName)
	maskFits = None
	if maskFile is not None:
		maskFits = fitsio.FITS(maskFile)
	statsPix = bokutil.stats_region(statreg)
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(4):
		ax = plt.subplot(2,2,i+1)
		extn = 'CCD%d'%(i+1)
		pix = fits[extn].read()[statsPix]
		if maskFits is not None:
			pix = np.ma.masked_array(pix,
			                         mask=maskFits[extn].read()[statsPix]>0)
		medVal,rms,pix = bokutil.array_stats(pix,method='median',
		                                     clip=True,rms=True,
		                                     retArray=True,**kwargs)
		meanVal = pix.mean()
		modalVal = 3*medVal-2*meanVal
		plt.hist(pix.flatten(),100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(medVal,color='r',lw=1.2)
		if showMean:
			plt.axvline(meanVal,color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,medVal,meanVal,
		         pix.std(),pix.min(),pix.max())

