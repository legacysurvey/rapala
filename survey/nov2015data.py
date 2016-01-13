#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip

from checkfields import srcor

default_zps = {'g':25.55,'r':25.23}

test_exptimes = np.array([100.]*6 + [25.,50.,100.,200.,400.])

def s82_ra_slice(s82,ra_range):
	return np.where((s82['ra']>ra_range[0]) & (s82['ra']<ra_range[1]))[0]

def load_stripe82_truth(ra_range=None):
	s82 = fits.getdata(os.path.join(os.environ['DESITARGTRUTH'],
	                                'stripe82-dr12-stars.fits.gz'))
	if ra_range is not None:
		s82 = s82[s82_ra_slice(s82,ra_range)]
	return s82

def has_coverage(truth,objs):
	return ( (truth['ra']>objs['ALPHA_J2000'].min()) &
	         (truth['ra']<objs['ALPHA_J2000'].max()) &
	         (truth['dec']>objs['DELTA_J2000'].min()) &
	         (truth['dec']<objs['DELTA_J2000'].max()) )

def match_stripe82():
	from collections import defaultdict
	bokdir = os.path.join(os.environ['BOKDATADIR'])
	s82 = load_stripe82_truth()
	cols = ['NUMBER','FLUX_APER','FLUXERR_APER','FLUX_AUTO','FLUXERR_AUTO',
	        'BACKGROUND','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000',
	        'ELONGATION','ELLIPTICITY','FWHM_IMAGE','FLAGS','CLASS_STAR',
	        'FLUX_RADIUS','FLUX_MAX','SNR_WIN','FLUX_PSF','FLUXERR_PSF']
	nAper = 3
	cdtypes = {'NUMBER':'i4','ALPHA_J2000':'f8','DELTA_J2000':'f8',
	           'FLAGS':'i8'}
	for filt in ['g','bokr']:
		s82_1 = s82_ra_slice(s82,(333.5,336.0))
		s82_2 = s82_ra_slice(s82,(351.0,353.5))
		fields1 = ['s82cal%s_ra334_%s' % (filt,n) for n in '123456abcde']
		fields2 = ['deep2%s_ra352_%s' % (filt,n) for n in '123456']
		for s82_ii,fields in zip([s82_1,s82_2],[fields1,fields2]):
			truth = s82[s82_ii]
			nTruth = len(s82_ii)
			nFields = len(fields)
			dtype = [('tIndex','i4'),('tMag','f4'),('tErr','f4'),
			         ('ccdNum','i4',(nFields,))]
			dtype.extend([(c,cdtypes.get(c,'f4'),(nFields,)) 
			                 for c in cols if 'APER' not in c])
			dtype.extend([(c,cdtypes.get(c,'f4'),(nFields,nAper)) 
			                 for c in cols if 'APER' in c])
			tab = np.zeros(nTruth,dtype=dtype)
			tab['tIndex'] = s82_ii
			tab['tMag'] = truth['psfMag_'+filt[-1]]
			tab['tErr'] = truth['psfMagErr_'+filt[-1]]
			tab['NUMBER'][:] = -1
			for fNum,field in enumerate(fields):
				imf = os.path.join(bokdir,filt,field+'.fits')
				catf = os.path.join(bokdir,filt,field+'.cat.fits')
				cat = fits.open(catf)
				for ccdNum in range(1,5):
					objs = cat[ccdNum].data
					ii = np.where(has_coverage(truth,objs))[0]
					m1,m2,d = srcor(objs['ALPHA_J2000'],objs['DELTA_J2000'],
					                truth['ra'][ii],truth['dec'][ii],2.0,
					                return_sep=True)
					print ' %20s[%d] matched %4d out of %4d (%.2f)' % \
					       (field,ccdNum,len(m1),len(ii),np.median(d))
					for col in cols:
						tab[col][ii[m2],fNum] = objs[col][m1]
					tab['ccdNum'][ii[m2],fNum] = ccdNum
				print
			tab = Table(tab)
			tab.write(fields[0][:-2]+'_merged.fits',overwrite=True)

def flux2mag(tab,band,which='PSF',zp=None):
	assert which in ['APER','AUTO','PSF']
	if zp is None:
		zp = default_zps[band] # XXX these are the nominal AM=1.3 vals
	k = 'FLUX_'+which
	mask = tab['NUMBER'] < 0
	if which=='APER':
		mask = np.dstack([mask]*tab[k].shape[-1])
	flux = np.ma.array(tab[k],mask=mask)
	exptimes = test_exptimes[:flux.shape[1]][np.newaxis,:]
	if which=='APER':
		exptimes = exptimes[:,:,np.newaxis]
	mag = zp - 2.5*np.ma.log10(flux/exptimes)
	return mag

def dump_zeropoints(tab,band):
	aperMag = flux2mag(tab,band,'APER')
	ref_star = (tab['tMag']>16.5) & (tab['tMag']<19.7)
	apNum = 2
	zp = default_zps[band]
	for j in range(tab['NUMBER'].shape[1]):
		print 'image %2d: ' % (j+1),
		for ccdNum in range(1,5):
			ii = np.where(~aperMag.mask[:,j,apNum] & ref_star &
			              (tab['ccdNum'][:,j]==ccdNum))[0]
			dm = aperMag[ii,j,apNum] - tab['tMag'][ii]
			dm = sigma_clip(dm)
			print '%3d %6.3f %6.3f  ' % \
			        (np.sum(~dm.mask),zp-dm.mean(),dm.std()),
		print

def dump_all_zeropoints():
	for fpfx,fsfx in [('s82cal','ra334'),('deep2','ra352')]:
		for filt in ['g','bokr']:
			fieldcatf = fpfx+filt+'_'+fsfx+'_merged.fits'
			fieldcat = Table.read(fieldcatf)
			print fieldcatf
			dump_zeropoints(fieldcat,filt[-1])
			print
		print
