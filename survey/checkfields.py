#!/usr/bin/env python

import os
import numpy as np
import fitsio

import bass
import bokextract

datadir = '/global/scratch2/sd/imcgreer/'
ndwfs_starfile = datadir+'ndwfs/starcat.fits'
cfhtlswide_starfile = datadir+'cfhtls/CFHTLSW3_starcat.fits'
cfhtlsdeep_starfile = datadir+'cfhtls/CFHTLSD3_starcat.fits'

def cfhtw3_tiles(observed=True):
	w3west,w3east = 15*(13.+50/60.), 15*(14+45./60)
	w3south,w3north = 50.7, 56.2
	return bass.region_tiles(w3west,w3east,w3south,w3north,observed=observed)

def ndwfs_tiles(observed=True):
	ndwest,ndeast = 15*14.37, 15*14.62
	ndsouth,ndnorth = 32.5, 36.1
	return bass.region_tiles(ndwest,ndeast,ndsouth,ndnorth,observed=observed)

def check_fields_list():
	files = [ t['utDate']+'/'+t['fileName']+'.fits.gz'
	                 for tiles in [cfhtw3_tiles(),ndwfs_tiles()] 
	                      for t in tiles ]
	with open('checkfields_tiles.txt','w') as f:
		f.write('\n'.join(sorted(files)))

def srcor(ra1,dec1,ra2,dec2,sep):
	from astropy.coordinates import SkyCoord,match_coordinates_sky
	from astropy import units as u
	c1 = SkyCoord(ra1,dec1,unit=(u.degree,u.degree))
	c2 = SkyCoord(ra2,dec2,unit=(u.degree,u.degree))
	idx,d2d,d3c = match_coordinates_sky(c1,c2)
	ii = np.where(d2d.arcsec < sep)[0]
	return ii,idx[ii]

def srcorXY(x1,y1,x2,y2,maxrad):
	sep = sqrt( (x1[:,np.newaxis]-x2[np.newaxis,:])**2 + 
	            (y1[:,np.newaxis]-y2[np.newaxis,:])**2 )
	ii = sep.argmin(axis=1)
	m1 = np.arange(len(x1))
	jj = np.where(sep[m1,ii] < maxrad)[0]
	return m1[jj],ii[jj]

def match_objects(objs,tiles):
	dtype = objs.dtype.descr + \
	           [('g_number','f4'),('g_ra','f8'),('g_dec','f8'),
	            ('g_autoMag','f4'),('g_autoMagErr','f4'),
	            ('g_autoFlux','f4'),('g_autoFluxErr','f4'),
	            ('g_elongation','f4'),('g_ellipticity','f4'),
	            ('g_flags','i4'),('g_fluxRad','f4'),
	            ('g_utDate','S8'),('g_expTime','f4'),
	            ('g_tileId','i4'),('g_ditherId','i4'),('g_ccdNum','i4')]
	skeys = ['NUMBER','ALPHA_J2000','DELTA_J2000','MAG_AUTO','MAGERR_AUTO',
	         'FLUX_AUTO','FLUXERR_AUTO','ELONGATION','ELLIPTICITY',
	         'FLAGS','FLUX_RADIUS']
	tkeys = ['utDate','expTime','tileId','ditherId']
	matches = []
	for ti,t in enumerate(tiles):
		print 'matching tile %d/%d' % (ti+1,len(tiles))
		for ccdNum in range(1,5):
			catpath = os.path.join(bass.rdxdir,t['utDate'],'ccdproc3',
			                       t['fileName']+'_ccd%d.cat.fits'%ccdNum)
			if not os.path.exists(catpath):
				print ' ... %s does not exist, skipping' % catpath
				continue
			cat = fitsio.read(catpath)
			ii = np.where( (objs['ra']>cat['ALPHA_J2000'].min()-1e-3) &
			               (objs['ra']<cat['ALPHA_J2000'].max()+1e-3) &
			               (objs['dec']>cat['DELTA_J2000'].min()-1e-3) &
			               (objs['dec']<cat['DELTA_J2000'].max()+1e-3) )[0]
			if len(ii)==0:
				continue
			m1,m2 = srcor(objs['ra'][ii],objs['dec'][ii],
			              cat['ALPHA_J2000'],cat['DELTA_J2000'],2.5)
			print '  ccd%d %d/%d' % (ccdNum,len(m1),len(ii)),
			matches.extend( [ tuple(objs[i]) +
			                  tuple([cat[k][j] for k in skeys]) +
			                  tuple([t[k] for k in tkeys]) + (ccdNum,)
			                     for i,j in zip(ii[m1],m2) ] )
		print ' total ',len(matches)
	matches = np.array(matches,dtype=dtype)
	print 'finished with ',matches.size
	return matches

##############################################################################
#                                                                            #
#                               NDWFS                                        #
#                                                                            #
##############################################################################

def select_ndwfs_stars():
	ndwfsdir = '/global/scratch2/sd/imcgreer/ndwfs/DR3/matchedFITS/'
	dtype = [('number','i4'),('autoMag','3f4'),('autoMagErr','3f4'),
	         ('ra','f8'),('dec','f8'),('rFWHM','f4'),('rClass','f4')]
	starcat = []
	rcols = ['NUMBER','MAG_AUTO','MAGERR_AUTO','ALPHA_J2000','DELTA_J2000',
	         'FWHM_IMAGE','CLASS_STAR']
	cols = ['MAG_AUTO','MAGERR_AUTO']
	for dec1 in range(32,36):
		catfn = lambda b: 'NDWFS_%s_%d_%d_cat_m.fits.gz' % (b,dec1,dec1+1)
		rfits = fitsio.FITS(ndwfsdir+catfn('R'))
		bfits = fitsio.FITS(ndwfsdir+catfn('Bw'))
		ifits = fitsio.FITS(ndwfsdir+catfn('I'))
		w = rfits[1].where('FWHM_IMAGE < 7 && MAG_AUTO < 24.7 && FLAGS == 0')
		print len(w)
		rcat = rfits[1].read(rows=w,columns=rcols)
		bcat = bfits[1].read(rows=w,columns=cols)
		icat = ifits[1].read(rows=w,columns=cols)
		stars = np.empty(len(w),dtype=dtype)
		stars['number'] = rcat['NUMBER']
		stars['ra'] = rcat['ALPHA_J2000']
		stars['dec'] = rcat['DELTA_J2000']
		stars['rFWHM'] = rcat['FWHM_IMAGE']
		stars['rClass'] = rcat['CLASS_STAR']
		for j,cat in enumerate([bcat,rcat,icat]):
			stars['autoMag'][:,j] = cat['MAG_AUTO']
			stars['autoMagErr'][:,j] = cat['MAGERR_AUTO']
		starcat.append(stars)
	starcat = np.concatenate(starcat)
	fitsio.write(ndwfs_starfile,starcat,clobber=True)

def match_ndwfs_stars(matchRad=2.5):
	stars = fitsio.read(ndwfs_starfile)
	tiles = ndwfs_tiles(observed=True)
	matches = match_objects(stars,tiles)
	fitsio.write('ndwfs_match.fits',matches,clobber=True)




##############################################################################
#                                                                            #
#                               CFHTLS                                       #
#                                                                            #
##############################################################################

def match_cfhtls_stars(matchRad=2.5,survey='wide'):
	if survey=='wide':
		stars = fitsio.read(cfhtlswide_starfile)
		tiles = cfhtw3_tiles(observed=True)
		fname = 'cfhtlswide'
	else:
		stars = fitsio.read(cfhtlsdeep_starfile)
		fname = 'cfhtlsdeep'
	matches = match_objects(stars,tiles)
	fitsio.write('%s_match.fits'%fname,matches,clobber=True)




##############################################################################
#                                                                            #
#                             fake sources                                   #
#                                                                            #
##############################################################################

from astropy.io import fits

def fake_sdss_stars_on_tile(stars,tile,
	                        nresample=100,magrange=(22.5,23.7),
	                        stampSize=25,margin=50,keepfakes=False):
	pixlo = lambda _x: _x-stampSize/2
	pixhi = lambda _x: _x-stampSize/2 + stampSize
	fakemags = np.zeros(nresample*4,dtype=np.float32)
	fakesnr = -np.ones_like(fakemags)
	for ccdNum in range(1,5):
		catpath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
		                       tile['fileName']+'_ccd%d.cat.fits'%ccdNum)
		if not os.path.exists(catpath):
			print ' ... %s does not exist, skipping' % catpath
			continue
		cat = fitsio.read(catpath)
		impath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
		                      tile['fileName']+'_ccd%d.fits'%ccdNum)
		fakeim = fits.open(impath)
		im = fakeim[0].data
		nY,nX = im.shape
		ii = np.where( (stars['ra']>cat['ALPHA_J2000'].min()+1e-3) &
		               (stars['ra']<cat['ALPHA_J2000'].max()-1e-3) &
		               (stars['dec']>cat['DELTA_J2000'].min()+1e-3) &
		               (stars['dec']<cat['DELTA_J2000'].max()-1e-3) )[0]
		if len(ii)==0:
			print 'no stars found on ccd #',ccdNum
			continue
		m1,m2 = srcor(stars['ra'][ii],stars['dec'][ii],
		              cat['ALPHA_J2000'],cat['DELTA_J2000'],2.5)
		jj = np.where(cat['FLAGS'][m2] == 0)[0]
		rindx = np.random.choice(len(jj),size=nresample,replace=True)
		fakemag = magrange[0] + \
		             (magrange[1]-magrange[0])*np.random.random(nresample)
		fscale = 10**(-0.4*(fakemag-stars['psfMag_g'][ii[m1[jj[rindx]]]]))
		fakex = np.random.randint(margin,nX-margin,nresample)
		fakey = np.random.randint(margin,nY-margin,nresample)
		for x,y,fx,fy,fscl in zip(cat['X_IMAGE'][m2[jj[rindx]]],
		                          cat['Y_IMAGE'][m2[jj[rindx]]],
		                          fakex,fakey,fscale):
			stamp = im[pixlo(y):pixhi(y),pixlo(x):pixhi(x)]
			im[pixlo(fy):pixhi(fy),pixlo(fx):pixhi(fx)] += fscl*stamp
		fakeimpath = impath.replace('.fits','_fake.fits')
		fakecatpath = fakeimpath.replace('.fits','.cat.fits')
		fakeim.writeto(fakeimpath,clobber=True)
		bokextract.sextract(fakeimpath)
		fakecat = fitsio.read(fakecatpath)
		q1,q2 = srcorXY(fakex,fakey,fakecat['X_IMAGE'],fakecat['Y_IMAGE'],3.0)
		snr = fakecat['FLUX_AUTO'][q2] / fakecat['FLUXERR_AUTO'][q2]
		fakemags[nresample*(ccdNum-1):nresample*ccdNum] = fakemag
		fakesnr[nresample*(ccdNum-1):nresample*ccdNum][q1] = snr
		if not keepfakes:
			os.unlink(fakeimpath)
			os.unlink(fakecatpath)
	return fakemags,fakesnr

def fake_ndwfs_stars(gmax=18.5,**kwargs):
	stars = fitsio.read('/global/scratch2/sd/imcgreer/ndwfs/sdss_bootes_gstars.fits')
	stars = stars[stars['psfMag_g']<gmax]
	tiles = ndwfs_tiles(observed=True)
	for ti,tile in enumerate(tiles):
		print 'faking stars in tile tile %d/%d' % (ti+1,len(tiles))
		fake_sdss_stars_on_tile(stars,tile,**kwargs)
		break

if __name__=='__main__':
	import sys
	if sys.argv[1]=='match_ndwfs':
		match_ndwfs_stars()
	elif sys.argv[1]=='match_cfhtlswide':
		print 'here'
		match_cfhtls_stars(survey='wide')

