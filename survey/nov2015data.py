#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.wcs import WCS

import matplotlib.pyplot as plt
from matplotlib import ticker

default_zps = {'g':25.55,'r':25.23}

test_exptimes = np.array([100.]*6 + [25.,50.,100.,200.,400.])

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

ampNums = [ [ 4,2,3,1 ] , [ 7,5,8,6 ], [ 10,12,9,11 ], [ 13,15,14,16] ]

def get_amp_index(tab):
	nx = 4096 // 2
	ny = 4032 // 2
	xi = (tab['X_IMAGE']/nx).astype(np.int32)
	yi = (tab['Y_IMAGE']/ny).astype(np.int32)
	ampIndex = 2*yi + xi
	return ampIndex

def calc_zeropoints(tab,band,apNum=2,zp=None,savedelta=False):
	if zp is None:
		zp = default_zps[band]
	aperMag = flux2mag(tab,band,'APER',zp=zp)
	ref_star = (tab['tMag']>16.5) & (tab['tMag']<19.7)
	nImages = tab['NUMBER'].shape[1]
	ampIndex = get_amp_index(tab)
	zpIm = np.zeros((nImages))
	zpCCD = np.zeros((nImages,4))
	zpAmp = np.zeros((nImages,4,4))
	for j in range(nImages):
		ii = np.where(~aperMag.mask[:,j,apNum] & ref_star)[0]
		dm = sigma_clip(aperMag[ii,j,apNum] - tab['tMag'][ii])
		if savedelta:
			zpIm[j] = dm.mean()
		else:
			zpIm[j] = zp - dm.mean()
		for ccdNum in range(1,5):
			ii = np.where(~aperMag.mask[:,j,apNum] & ref_star &
			              (tab['ccdNum'][:,j]==ccdNum))[0]
			dm = sigma_clip(aperMag[ii,j,apNum] - tab['tMag'][ii])
			if savedelta:
				zpCCD[j,ccdNum-1] = dm.mean()
			else:
				zpCCD[j,ccdNum-1] = zp - dm.mean()
			for ai in range(4):
				ii = np.where(~aperMag.mask[:,j,apNum] & ref_star &
				              (tab['ccdNum'][:,j]==ccdNum) & 
				              (ampIndex[:,j]==ai))[0]
				dm = sigma_clip(aperMag[ii,j,apNum] - tab['tMag'][ii])
				if savedelta:
					zpAmp[j,ccdNum-1,ai] = dm.mean()
				else:
					zpAmp[j,ccdNum-1,ai] = zp - dm.mean()
	return zpIm,zpCCD,zpAmp

def dump_zeropoints(tab,band,byamp=False,**kwargs):
	zpim,zpccd,zpamp = calc_zeropoints(tab,band,**kwargs)
	for j in range(len(zpim)):
		print 'image %2d: ' % (j+1),
		for ccdNum in range(1,5):
			print '%6.3f ' % zpccd[j,ccdNum-1],
			if not byamp:
				continue
			for ai in range(4):
				print '  IM%d  %6.3f  ' % \
				   (ampNums[ccdNum-1][ai],zpamp[j,ccdNum-1,ai]),
			print
		print '   %6.3f' % zpim[j]

def delta_zps_byamp(zpFit):
	for ccdi in range(4):
		ccdmean = np.mean(zpFit[:,ccdi,:])
		for ampj in range(4):
			ampmean = np.mean(zpFit[:,ccdi,ampj])
			print 'IM%d %8.4f' % (ampNums[ccdi][ampj],ccdmean-ampmean)

def dump_all_zeropoints():
	for fpfx,fsfx in [('s82cal','ra334'),('deep2','ra352')]:
		for filt in ['g','bokr']:
			fieldcatf = fpfx+filt+'_'+fsfx+'_merged.fits'
			fieldcat = Table.read(fieldcatf)
			print fieldcatf
			dump_zeropoints(fieldcat,filt[-1])
			print
		print

def check_scatter(tab,band):
	from scipy.stats import scoreatpercentile
	psfMag = flux2mag(tab,band,'PSF')
	for j in range(tab['NUMBER'].shape[1]):
		print 'image %2d: ' % (j+1)
		for mag1 in np.arange(17,21.1,0.5):
			is_mag = (tab['tMag']>(mag1-0.25)) & (tab['tMag']<(mag1+0.25))
			ii = np.where(~psfMag.mask[:,j] & is_mag)[0]
			dm = psfMag[ii,j] - tab['tMag'][ii]
			dist = scoreatpercentile(np.abs(dm-np.median(dm)),
			                         [25,50,75,90,95])
			print '  ',mag1,dist

def focalplanevar(tab,band,nbin=4,doplot=False,vr=0.015,shownum=False,
                  zpcorr=None):
	apNum = 2
	if zpcorr is None:
		zp = None
	else:
		zpim,zpccd,zpamp = calc_zeropoints(tab,band,apNum=apNum)
		if zpcorr=='image':
			zp = zpim[:,np.newaxis]
		elif zpcorr=='ccd':
			# add zp=0 for ccd=0 (empty entries)
			zp = np.hstack([np.zeros((zpccd.shape[0],1)),zpccd])
			zp = np.choose(tab['ccdNum'],zp.transpose())[:,:,np.newaxis]
		elif zpcorr=='amp':
			ampIndex = get_amp_index(tab)
			ai = np.clip(4*(tab['ccdNum']-1) + ampIndex + 1, 0, 16)
			zp = np.hstack([np.zeros((zpamp.shape[0],1)),zpamp.reshape(-1,16)])
			zp = np.choose(ai,zp.transpose())[:,:,np.newaxis]
	aperMag = flux2mag(tab,band,'APER',zp=zp)[:,:,apNum]
	refMag = tab['tMag']
	meanMag = sigma_clip(aperMag,axis=1).mean(axis=1)
	deltaMag = aperMag - meanMag[:,np.newaxis]
	deltaMag[tab['tMag']>19,:].masked = True
	nx = 4096 // nbin
	ny = 4032 // nbin
	xi = (tab['X_IMAGE']/nx).astype(np.int32)
	yi = (tab['Y_IMAGE']/ny).astype(np.int32)
	fpIm = np.zeros((4,nbin,nbin))
	for ccdi in range(4):
		for i in range(nbin):
			for j in range(nbin):
				ii = np.where((tab['ccdNum']==ccdi+1)&(yi==i)&(xi==j))
				dm = sigma_clip(deltaMag[ii])
				if shownum:
					fpIm[ccdi,i,j] = (~dm.mask).sum()
				else:
					fpIm[ccdi,i,j] = np.ma.mean(dm)
	if doplot:
		fig = plt.figure(figsize=(6,6.15))
		plt.subplots_adjust(0.04,0.035,0.96,0.88,0.25,0.12)
		for pnum,ccdi in enumerate([0,2,1,3],start=1):
			ax = plt.subplot(2,2,pnum)
			im = fpIm[ccdi]
			if ccdi <= 1:
				im = im[:,::-1]
			if ccdi % 2 == 1:
				im = im[::-1,:]
			if shownum:
				_im = ax.imshow(im,origin='lower',interpolation='nearest',
				                cmap=plt.cm.hot_r)
			else:
				_im = ax.imshow(im,origin='lower',interpolation='nearest',
				                vmin=-vr,vmax=vr,cmap=plt.cm.RdBu)
			plt.title('CCD%d'%(ccdi+1))
			if pnum==1:
				cbax = fig.add_axes([0.1,0.98,0.8,0.015])
				cb = fig.colorbar(_im,cax=cbax,orientation='horizontal')
				if not shownum:
					cb.locator = ticker.MultipleLocator(0.005)
				cb.update_ticks()
			ax.xaxis.set_ticks([])
			ax.yaxis.set_ticks([])
	return fpIm

def stripe82zps():
	s82g = Table.read('s82calg_ra334_merged.fits')
	s82r = Table.read('s82calbokr_ra334_merged.fits')
	plt.ioff()
	for filt,tab in [('g',s82g),('r',s82r)]:
		for zpcorr in ['image','ccd','amp']:
			fpim = focalplanevar(tab,filt,doplot=True,nbin=8,vr=0.015,
			                     zpcorr=zpcorr)
			plt.savefig('s82zp_%s_%s.png' % (filt,s))
	plt.ion()

def stripe82_depth(imageFile,aperRad=7.0):
	from bokpipe.bokphot import aper_phot_image
	_minmax = lambda a,b: (a,b) if a<=b else (b,a)
	aperRad /= 0.455
	s82all = load_stripe82_truth(ra_range=(332,336))
	ph = aper_phot_image(imageFile,s82all['ra'],s82all['dec'],[aperRad])
	ph['refMag'] = s82all['psfMag_g'][ph['objId']]
	ph['ra'] = s82all['ra'][ph['objId']]
	ph['dec'] = s82all['dec'][ph['objId']]
	return ph

