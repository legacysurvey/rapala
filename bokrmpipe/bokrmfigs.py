#!/usr/bin/env python

import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy.table import Table
from astropy.io import fits

def zp_strip_chart(kind='aper',byFrame=False):
	plt.figure(figsize=(14,4))
	plt.subplots_adjust(0.04,0.10,0.98,0.98)
	for pNum,filt in enumerate('g',start=1):
		ax = plt.subplot(1,1,pNum)
		zpDat = Table.read('zeropoints_%s.fits'%filt)
		utds,ii = np.unique(zpDat['utDate'],return_index=True)
		xval = zpDat['frameNum'] if byFrame else np.arange(len(zpDat))
		plt.plot(xval,zpDat[kind+'Zp'])
		for i,utd in zip(ii,utds):
			plt.axvline(xval[i],c='gray')
			plt.text(xval[i]+1,25.55,utd,
			         ha='left',rotation='vertical',size=9)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
		ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
		ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
	plt.xlim(xval[0]-20,xval[-1]+20)
	plt.ylim(25.2,26.15)

def plot_lightcurve(targetNum,targetSource='RM',shownightly=False):
	pfx = 'bokrm_sdss'
	gcat = Table.read('lightcurves_bokrm_g.fits')
	icat = None
	refCat = fits.getdata(os.environ['BOK90PRIMEDIR']+'/../data/sdss.fits',1)
	aperNum = 3
	ymin,ymax = 1e9,0
	plt.figure(figsize=(12,5))
	plt.subplots_adjust(0.05,0.05,0.97,0.94)
	ax,pnum = None,1
	gcat = gcat.group_by('objId')
	j = np.where(gcat.groups.keys['objId']==targetNum)[0][0]
	glc = gcat.groups[j]
	ilc = None #icat.groups[j]
	for band,lc,clr in zip('gi',[glc,ilc],['g','r']):
		ax = plt.subplot(2,1,pnum,sharex=ax)
		jj = np.where(lc['aperMagErr'][:,aperNum] > 0)[0]
		if len(jj)==0:
			pnum += 1
			continue
		plt.errorbar(lc['mjd'][jj],
		             lc['aperMag'][jj,aperNum],
		             lc['aperMagErr'][jj,aperNum],
		             fmt='s',mfc='none',ecolor=clr,ms=2)
		if shownightly:
			nightly = Table.read('nightly_lcs_bokrm_g.fits')
			nightly = nightly.group_by('objId')
			_j = np.where(nightly.groups.keys['objId']==targetNum)[0][0]
			_lc = nightly.groups[j]
			plt.errorbar(_lc['mean_mjd'][:],
			             _lc['aperMag'][:,aperNum],
			             _lc['aperMagErr'][:,aperNum],
			             fmt='s',mfc='none',mew=2,ecolor=clr,ms=8)
		ymin = np.percentile(lc['aperMag'][jj,aperNum],10)
		ymax = np.percentile(lc['aperMag'][jj,aperNum],90)
		dy = max(0.05,5*np.median(lc['aperMagErr'][jj,aperNum]))
		plt.ylim(ymin-dy,ymax+dy)
		bi = 1 if clr=='g' else 3
		if False:#targetSource=='RM':
			plt.axhline(target['PSFMAG'][targetNum,bi],color=clr)
			plt.scatter(target['DR_MJD'][targetNum],
			            target['PSFMAG'][targetNum,bi],
			            c=clr)
			if clr=='g':
				plt.title('rm%03d, z=%.2f' % 
				          (targetNum,target['ZFINAL'][targetNum]))
		pnum += 1
		break
	# XXX change lims when have i band
	plt.xlim(min(glc['mjd'].min(),np.inf)-5,max(glc['mjd'].max(),0)+5)

