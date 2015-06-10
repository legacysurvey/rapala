#!/usr/bin/env python

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
	ff = fitsio.FITS('/global/scratch2/sd/imcgreer/ndwfs/starcat.fits','rw')
	ff.write(starcat,clobber=True)

