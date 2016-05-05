#!/usr/bin/env python

import os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,hstack,vstack
from astropy.stats import sigma_clip
import healpy

class NoCalibrationStars(Exception):
	pass

def srcor(ra1,dec1,ra2,dec2,sep,return_sep=False):
	from astropy.coordinates import SkyCoord,match_coordinates_sky
	from astropy import units as u
	c1 = SkyCoord(ra1,dec1,unit=(u.degree,u.degree))
	c2 = SkyCoord(ra2,dec2,unit=(u.degree,u.degree))
	idx,d2d,d3c = match_coordinates_sky(c1,c2)
	ii = np.where(d2d.arcsec < sep)[0]
	if return_sep:
		return ii,idx[ii],d2d.arcsec[ii]
	else:
		return ii,idx[ii]

def read_ps1cat(ra,dec):
	#
	path = '/global/project/projectdirs/cosmo/work/ps1/cats/chunks-qz-star-v2/'
	#
	phi = np.radians(ra)
	theta = np.radians(90-dec)
	pixels = healpy.pixelfunc.ang2pix(32,theta,phi)
	pixels = np.unique(pixels)
	#
	cat = []
	for pix in pixels:
		fname = path+'ps1-%05d.fits'%pix
		if os.path.exists(fname):
			cat.append(fits.getdata(fname,1))
	if not cat:
		raise NoCalibrationStars
	else:
		return np.concatenate(cat)

def get_ps1_stars(ra,dec):
	ps1cat = read_ps1cat(ra,dec)
	gicolor = ps1cat['MEDIAN'][:,0] - ps1cat['MEDIAN'][:,2]
	delta = 20./3600
	ii = np.where( (gicolor>0.4) & (gicolor<2.7) &
	               (ps1cat['NMAG_OK'][:,0] >= 1) &
	               (ps1cat['NMAG_OK'][:,2] >= 1) &
	               (ps1cat['RA']<ra.max()+delta) &
	               (ps1cat['RA']>ra.min()-delta) &
	               (ps1cat['DEC']<dec.max()+delta) &
	               (ps1cat['DEC']>dec.min()-delta) )[0]
	return ps1cat[ii]

def match_ps1(catf,stars=True,isldac=False,radius=10.):
	cats = []
	fitscat = fits.open(catf)
	for ccdNum in range(1,5):
		if isldac:
			cat = fitscat[2*ccdNum].data
		else:
			cat = fitscat[ccdNum].data
		try:
			if stars:
				ps1objs = get_ps1_stars(cat['ALPHA_J2000'],cat['DELTA_J2000'])
			else:
				ps1objs = read_ps1cat(cat['ALPHA_J2000'],cat['DELTA_J2000'])
		except NoCalibrationStars:
			continue
		if len(ps1objs)<10:
			continue
		m1,m2,d = srcor(ps1objs['RA'],ps1objs['DEC'],
		                cat['ALPHA_J2000'],cat['DELTA_J2000'],
		                radius,return_sep=True)
		#print catf,ccdNum,len(cat),len(ps1objs),len(m1)
		t2 = Table(cat[m2])
		t2['ccdNum'] = ccdNum
		t2['separation'] = d
		cats.append(hstack([Table(ps1objs[m1]),t2]))
	if not cats:
		return None
	return vstack(cats)

