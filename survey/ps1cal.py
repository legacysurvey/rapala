#!/usr/bin/env python

import os
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table,hstack,vstack
from astropy.stats import sigma_clip
import healpy

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
	return np.concatenate(cat)

def get_ps1_stars(ra,dec):
	ps1cat = read_ps1cat(ra,dec)
	gicolor = ps1cat['MEDIAN'][:,0] - ps1cat['MEDIAN'][:,2]
	delta = 20./3660
	ii = np.where( (gicolor>0.4) & (gicolor<2.7) &
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
		if stars:
			ps1objs = get_ps1_stars(cat['ALPHA_J2000'],cat['DELTA_J2000'])
		else:
			ps1objs = read_ps1cat(cat['ALPHA_J2000'],cat['DELTA_J2000'])
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
	if len(cats)==0:
		return None
	return vstack(cats)

def match_nov15data():
	#bokdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',
	#                      'bokpipe_v0.2','nov15data')
	bokdir = '/global/homes/i/imcgreer/bok/reduced/nov15data_ian/'
	for filt in 'gr':
		fields = ['s82cal%s_ra334_%s' % (filt,n) for n in '123456abcde']
		fields += ['deep2%s_ra352_%s' % (filt,n) for n in '123456']
		for field in fields:
			imf = os.path.join(bokdir,filt,field+'.fits')
			catf = imf.replace('.fits','.cat.fits')
			cat = match_ps1(catf)
			print field,len(cat)
			cat.write(os.path.join('ps1qz',field+'_ps1qz.fits'))

def nov15_ps1_zps(apNum=2,doplots=False):
	import nov2015data
	fluxk = 'FLUX_APER'
	bokdir = '/global/homes/i/imcgreer/bok/reduced/nov15data_ian/'
	d = defaultdict(list)
	for j,band in enumerate('gr'):
		fields = ['s82cal%s_ra334_%s' % (band,n) for n in '123456abcde']
		fields += ['deep2%s_ra352_%s' % (band,n) for n in '123456']
		exptimes = np.concatenate([nov2015data.test_exptimes,
		                           nov2015data.test_exptimes[:-5]])
		zp = nov2015data.bok_zpt0[{'g':'g','r':'bokr'}[band]]
		for field,texp in zip(fields,exptimes):
			cat = Table.read(os.path.join('ps1qz',field+'_ps1qz.fits'))
			good = ( (cat['MEDIAN'][:,j] < 25) & (cat['MEDIAN'][:,j]>16) &
			         (cat['FLAGS'] == 0) )
			aperMag = zp - 2.5*np.log10(cat[fluxk]/texp)
			ii = np.where(good)[0]
			dm = sigma_clip(aperMag[ii,apNum] - cat['MEDIAN'][ii,j])
			zpIm = zp - dm.mean()
			zpCCD = np.zeros(4)
			zpAmp = np.zeros((4,4))
			ampIndex = nov2015data.get_amp_index(cat['X_IMAGE'],cat['Y_IMAGE'])
			for ccdNum in range(1,5):
				ii = np.where(good & (cat['ccdNum']==ccdNum))[0]
				dm = sigma_clip(aperMag[ii,apNum] - cat['MEDIAN'][ii,j])
				zpCCD[ccdNum-1] = zp - dm.mean()
				for ai in range(4):
					ii = np.where(good &
					              (cat['ccdNum']==ccdNum) & 
					              (ampIndex==ai))[0]
					dm = sigma_clip(aperMag[ii,apNum] - cat['MEDIAN'][ii,j])
					zpAmp[ccdNum-1,ai] = zp - dm.mean()
			d['fileName'].append(field)
			d['zpIm'].append(zpIm)
			d['zpCCD'].append(zpCCD)
			d['zpAmp'].append(zpAmp)
	tab = Table(d)
	tab.write('nov2015_ps1_zps.fits',overwrite=True)

def ps1colorterms(doplots=False):
	import nov2015data
	fluxk = 'FLUX_APER'
	bokdir = '/global/homes/i/imcgreer/bok/reduced/nov15data_ian/'
	zptab = Table.read('nov2015_ps1_zps.fits')
	zptab.add_index('fileName')
	apNum=2 # XXX
	if doplots:
		plt.figure(figsize=(12,5))
		plt.subplots_adjust(0.10,0.14,0.97,0.92,0.31,0.2)
	for j,band in enumerate('gr'):
		if doplots:
			plt.subplot(1,2,j+1)
		fields = ['s82cal%s_ra334_%s' % (band,n) for n in '123456abcde']
		fields += ['deep2%s_ra352_%s' % (band,n) for n in '123456']
		exptimes = np.concatenate([nov2015data.test_exptimes,
		                           nov2015data.test_exptimes[:-5]])
		gicolor = []
		dmag = []
		for field,texp in zip(fields,exptimes):
			cat = Table.read(os.path.join('ps1qz',field+'_ps1qz.fits'))
			good = ( (cat['MEDIAN'][:,j] < 25) & (cat['MEDIAN'][:,j]>16) &
			         (cat['FLAGS'] == 0) )
			ampIndex = nov2015data.get_amp_index(cat['X_IMAGE'],cat['Y_IMAGE'])
			zp = zptab.loc[field]['zpAmp'][cat['ccdNum']-1,ampIndex]
			aperMag = zp - 2.5*np.log10(cat[fluxk][:,apNum]/texp)
			gicolor.append(np.diff(cat['MEDIAN'][:,[2,0]],axis=1).squeeze())
			dmag.append(aperMag - cat['MEDIAN'][:,j])
		gicolor = np.concatenate(gicolor)
		dmag = np.concatenate(dmag)
		if doplots:
			ii = np.where(np.abs(dmag) < 0.2)[0]
			plt.hexbin(gicolor[ii],dmag[ii],cmap=plt.get_cmap('gray_r'),
			           bins='log')
			plt.xlabel('PS1 g-i')
			plt.ylabel('%s(Bok) - %s(PS1)' % (band,band))
		order = 3
		mask = np.abs(dmag) > 0.25
		for iternum in range(3):
			order = iternum+1
			_dmag = np.ma.array(dmag,mask=mask)
			fit = np.ma.polyfit(gicolor,_dmag,order)
			magoff = sigma_clip(dmag-np.polyval(fit,gicolor))#,iters=1)
			mask = magoff.mask
			if doplots and iternum==2:
				xx = np.linspace(0.4,2.7,100)
				plt.plot(xx,np.polyval(fit,xx),c='gbr'[iternum])
		if doplots:
			polystr = ' '.join(['%+.5f*gi^%d'%(c,order-d) 
			                      for d,c in enumerate(fit)])
			plt.title(('$%s(Bok) - %s(PS1) = '%(band,band))+polystr+'$',
			          size=11)
		np.savetxt('bok2ps1_%s_coeff.txt'%band,fit)
	if doplots:
		plt.savefig('bok2ps1.png')

