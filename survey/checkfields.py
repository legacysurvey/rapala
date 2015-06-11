#!/usr/bin/env python

import os
import numpy as np
import fitsio

import bass

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

ndwfs_starfile = '/global/scratch2/sd/imcgreer/ndwfs/starcat.fits'

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

def srcor(ra1,dec1,ra2,dec2,sep):
	from astropy.coordinates import SkyCoord,match_coordinates_sky
	from astropy import units as u
	c1 = SkyCoord(ra1,dec1,unit=(u.degree,u.degree))
	c2 = SkyCoord(ra2,dec2,unit=(u.degree,u.degree))
	idx,d2d,d3c = match_coordinates_sky(c1,c2)
	ii = np.where(d2d.arcsec < sep)[0]
	return ii,idx[ii]

def match_ndwfs_stars(matchRad=2.5):
	stars = fitsio.read(ndwfs_starfile)
	tiles = ndwfs_tiles(observed=True)
	dtype = stars.dtype.descr + \
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
			ii = np.where( (stars['ra']>cat['ALPHA_J2000'].min()-1e-3) &
			               (stars['ra']<cat['ALPHA_J2000'].max()+1e-3) &
			               (stars['dec']>cat['DELTA_J2000'].min()-1e-3) &
			               (stars['dec']<cat['DELTA_J2000'].max()+1e-3) )[0]
			if len(ii)==0:
				continue
			m1,m2 = srcor(stars['ra'][ii],stars['dec'][ii],
			              cat['ALPHA_J2000'],cat['DELTA_J2000'],2.5)
			print '  ccd%d %d/%d' % (ccdNum,len(m1),len(ii)),
			matches.extend( [ tuple(stars[i]) +
			                  tuple([cat[k][j] for k in skeys]) +
			                  tuple([t[k] for k in tkeys]) + (ccdNum,)
			                     for i,j in zip(ii[m1],m2) ] )
		print ' total ',len(matches)
	matches = np.array(matches,dtype=dtype)
	print 'finished with ',matches.size
	fitsio.write('ndwfs_match.fits',matches,clobber=True)

if __name__=='__main__':
	import sys
	if sys.argv[1]=='match_ndwfs':
		match_ndwfs_stars()

