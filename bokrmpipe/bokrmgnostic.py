#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astropy import units as u

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

def srcor(ra1,dec1,ra2,dec2,sep):
	c1 = SkyCoord(ra1,dec1,unit=(u.degree,u.degree))
	c2 = SkyCoord(ra2,dec2,unit=(u.degree,u.degree))
	idx,d2d,d3c = match_coordinates_sky(c1,c2)
	ii = np.where(d2d.arcsec < sep)[0]
	return ii,idx[ii],d2d.arcsec[ii]

def check_img_astrom(imgFile,refCat,catFile=None,mlim=19.5,band='g'):
	imFits = fits.open(imgFile)
	if catFile is None:
		catFile = imgFile.replace('.fits','.cat.fits')
	catFits = fits.open(catFile)
	try:
		ahead = bokastrom.read_headers(imgFile.replace('.fits','.ahead'))
	except:
		ahead = None
	rv = []
	for ccd in range(1,5):
		ccdCat = catFits[ccd].data
		hdr = imFits[ccd].header
		if ahead is not None:
			hdr.update(ahead[ccd-1].items())
		w = WCS(hdr)
		foot = w.calc_footprint()
		ras = sorted(foot[:,0])
		decs = sorted(foot[:,1])
		ii = np.where((refCat['ra']>ras[1])&(refCat['ra']<ras[2]) &
		              (refCat['dec']>decs[1])&(refCat['dec']<decs[2]) &
		              (refCat[band]<mlim))[0]
		m1,m2,sep = srcor(ccdCat['ALPHA_J2000'],ccdCat['DELTA_J2000'],
		                  refCat['ra'][ii],refCat['dec'][ii],5.0)
		rv.append(dict(N=len(ii),nMatch=len(ii),
		               ra=ccdCat['ALPHA_J2000'][m1],
		               dec=ccdCat['DELTA_J2000'][m1],
		               raRef=refCat['ra'][ii[m2]],
		               decRef=refCat['dec'][ii[m2]],
		               sep=sep))
	return rv

def check_astrom(files):
	sdss = fits.getdata(os.environ['BOK90PRIMEDIR']+'/../data/sdss.fits',1)
	for f in files:
		m = check_img_astrom(f,sdss)
		print os.path.basename(f),
		for c in m:
			print '%7.3f' % np.median(c['sep']),
		print

