#!/usr/bin/env python

import os
import itertools
import numpy as np
from astropy.io import fits
from astropy.table import Table,join,vstack
from astropy.stats import sigma_clip
from astropy.wcs import WCS
import fitsio

import matplotlib.pyplot as plt
from matplotlib import ticker

from bokpipe.bokphot import aper_phot_image

from bass import ampNums,get_amp_index
import ps1cal

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

def load_obsdb():
	return Table.read('../basschute/config/nov2015_mod.fits')

def stripe82_catalogs(bokdir,ccdsfile='bass-ccds-idmnov2015.fits'):
	from boketc import k_ext
	alltabs = []
	mrows = []
	fieldNum = 1
	ccdlist = Table.read(ccdsfile)
	# start with the Stripe 82 truth table as a base catalog
	s82truth = load_stripe82_truth()
	# iterate over the two S82 fields observed (Deep2F3 and a random spot)
	for field,ra_range in zip(['s82cal%s_ra334_%s','deep2%s_ra352_%s'],
	                          [(333.5,336.0),(351.0,353.5)]):
		ii = s82_ra_slice(s82truth,ra_range)
		s82 = s82truth[ii]
		# initialize the table with Stripe 82 fields and do some renaming
		# for convenience
		t = Table(s82truth[ii],masked=True)
		for cn in t.colnames:
			cn_new = cn.replace('psfmagerr','err')
			cn_new = cn_new.replace('psfmag','mag')
			t.rename_column(cn,'sdss_'+cn_new)
		t['sdss_idx'] = ii
		# match the Bok catalogs
		for k,dtyp in [('x','f4'),('y','f4'),('ra','f8'),('dec','f8'),
		               ('idx','i4'),('mag','f4'),('err','f4'),
		               ('field','i4'),('ccd','i4'),('amp','i4'),
		               ('kx','f4'),('apcorr','f4')]:
			for b in 'gr':
				t['bok_'+k+'_'+b] = np.zeros((len(ii),11),dtype=dtyp)
		for b in 'gr':
			t['bok_mag_'+b].mask = True
		for n,b in itertools.product(range(11),'gr'):
			sfx = '123456abcde'[n]
			fieldName = field % (b,sfx)
			catf = os.path.join(bokdir,fieldName+'.cat.fits')
			try:
				cat = fits.open(catf)
			except IOError:
				continue
			for ccdNum in range(1,5):
				objs = cat[ccdNum].data
				ii = np.where(has_coverage(s82,objs))[0]
				m1,m2,d = srcor(objs['ALPHA_J2000'],objs['DELTA_J2000'],
				                s82['ra'][ii],s82['dec'][ii],
				                2.0,return_sep=True)
				print ' %20s[%d] matched %4d out of %4d (%.2f)' % \
				       (fieldName,ccdNum,len(m1),len(ii),np.median(d))
				t['bok_idx_'+b][ii[m2],n] = objs['NUMBER'][m1]
				t['bok_field_'+b][ii[m2],n] = fieldNum
				t['bok_ccd_'+b][ii[m2],n] = ccdNum
				ampIndex = get_amp_index(objs['X_IMAGE'][m1],
				                         objs['Y_IMAGE'][m1])
				t['bok_amp_'+b][ii[m2],n] = 4*(ccdNum-1) + ampIndex + 1
				for k1,k2 in [('x','X_IMAGE'),('y','Y_IMAGE'),
				              ('ra','ALPHA_J2000'),('dec','DELTA_J2000')]:
					t['bok_'+k1+'_'+b][ii[m2],n] = objs[k2][m1]
				flux = np.ma.array(objs['FLUX_APER'][m1,1],
				                   mask=((objs['FLUX_APER'][m1,1]<=0) |
				                         (objs['FLAGS'][m1]>0)))
				fluxerr = objs['FLUXERR_APER'][m1,1]
				mag = -2.5*np.ma.log10(flux/test_exptimes[n])
				t['bok_mag_'+b][ii[m2],n] = mag
				t['bok_mag_'+b][ii[m2],n].mask = mag.mask
				t['bok_err_'+b][ii[m2],n] = 1.0856*np.ma.divide(fluxerr,flux)
				# meta-data about the image
				totflux = objs['FLUX_APER'][m1,-1]
				apcorr = objs['FLUX_APER'][m1] / totflux[:,np.newaxis]
				psfcorr = objs['FLUX_PSF'][m1] / totflux
				fwhm = np.median(objs['FWHM_IMAGE'][m1])
				apcorr = sigma_clip(apcorr[:,1]).mean()
				ccdent = np.where((ccdlist['image_filename'] ==
			                                 fieldName+'.fits') &
				                        (ccdlist['ccdnum']==ccdNum))[0][0]
				airmass = ccdlist['airmass'][ccdent]
				mrows.append((fieldName,fieldNum,b,ccdNum,fwhm,apcorr,
				              airmass,ccdlist['avsky'][ccdent],
				              ccdlist['arawgain'][ccdent]))
				# copy corrections into photo tab
				t['bok_kx_'+b][ii[m2],n] = k_ext[b]*(airmass-1)
				t['bok_apcorr_'+b][ii[m2],n] = 2.5*np.log10(apcorr)
			if b=='r':
				fieldNum += 1
		# eliminate entries without Bok coverage in both bands
		ii = np.where(~(t['bok_mag_g'].mask.all(axis=1) &
		                t['bok_mag_r'].mask.all(axis=1)))[0]
		t = t[ii]
		# now bring in PS1
		for k,dtyp in [('ra','f8'),('dec','f8'),('idx','i4')]:
			t['ps1_'+k] = np.zeros(len(t),dtype=dtyp)
		for j,b in enumerate('grizy'):
			t['ps1_mag_'+b] = np.zeros(len(t),dtype='f4')
		ps1objs = ps1cal.get_ps1_stars(t['sdss_ra'],t['sdss_dec'])
		m1,m2 = srcor(ps1objs['RA'],ps1objs['DEC'],
		              t['sdss_ra'],t['sdss_dec'],2.0)
		for k1,k2 in [('ra','RA'),('dec','DEC'),('idx','OBJ_ID')]:
			t['ps1_'+k1][m2] = ps1objs[k2][m1]
		for j,b in enumerate('grizy'):
			t['ps1_mag_'+b][m2] = ps1objs['MEDIAN'][m1,j]
		# 
		alltabs.append(t)
	# combine tables & calculate residual per-amplifier offsets
	tab = vstack(alltabs)
	for b in 'gr':
		mag = tab['bok_mag_'+b] - tab['bok_kx_'+b] + tab['bok_apcorr_'+b]
		mivar = np.ma.power(tab['bok_err_'+b],-2)
		mean_mag = np.ma.average(sigma_clip(mag,axis=1),weights=mivar,axis=1)
		tab['bok_meanmag_'+b] = mean_mag
		tab['bok_ampcorr_'+b] = np.zeros_like(tab['bok_mag_'+b])
		print '\n%s band amp correction (millimag): ' % b
		for amp in range(1,17):
			amp_mags = np.ma.array(mag,mask=(mag.mask|tab['bok_amp_'+b]!=amp))
			amp_dmag = sigma_clip(mean_mag[:,np.newaxis] - amp_mags)
			ii = np.where(tab['bok_amp_'+b] == amp)
			tab['bok_ampcorr_'+b][ii] = amp_dmag.mean()
			tab.meta['AMPC%s%02d'%(b.upper(),amp)] = amp_dmag.mean()
			print '%4d +/- %2d ' % (1e3*amp_dmag.mean(),1e3*amp_dmag.std()),
			if (amp%4)==0: print
		print
		mag += tab['bok_ampcorr_'+b]
		tab['bok_magcorr_'+b] = mag 
		mean_mag = np.ma.average(sigma_clip(mag,axis=1),weights=mivar,axis=1)
		tab['bok_meanmagcorr_'+b] = mean_mag 
	# something weird here
	for k in ['mag','magcorr','ampcorr']:
		for b in 'gr':
			tab['bok_'+k+'_'+b].fill_value = np.nan
	tab.write('stripe82bok_nov2015stars.fits',overwrite=True)
	#
	metaTab = Table(rows=mrows,
	                names=('field','fieldNum','filter','ccdNum',
	                       'fwhm','apCorr','airmass','skyADU','gain'))
	metaTab.write('stripe82bok_nov2015sum.fits',overwrite=True)

def load_nov2015_data():
	t = Table.read('stripe82bok_nov2015stars.fits')
	# multi-dim Table columns written to FITS lose their masks...
	for k in ['mag','magcorr','meanmag','meanmagcorr']:
		for b in 'gr':
			t['bok_%s_%s'%(k,b)].mask = np.isnan(t['bok_%s_%s'%(k,b)])
	for b in 'grizy':
		t['ps1_mag_'+b].mask = t['ps1_mag_'+b]==0
	return t

def calc_color_terms(doplots=False,savefit=False):
	t = load_nov2015_data()
	n = t['bok_mag_g'].shape[1]
##	refclr = {'sdss':np.tile(t['sdss_mag_g']-t['sdss_mag_z'],(n,1)).transpose(),
##	           'ps1':np.tile(t['ps1_mag_g']-t['ps1_mag_i'],(n,1)).transpose()}
	refclr = {'sdss':t['sdss_mag_g']-t['sdss_mag_i'],
	           'ps1':t['ps1_mag_g']-t['ps1_mag_i']}
	refs = ['sdss','ps1']
	for b in 'gr':
		refclr_all = {'sdss':[],'ps1':[]}
		dmag_all = {'sdss':[],'ps1':[]}
		mag = t['bok_meanmagcorr_'+b]
		for ref in refs:
			refmag = t[ref+'_mag_'+b]
			dmag = sigma_clip(mag-refmag)
			zp0 = dmag.mean()
			bokmag = mag - zp0
			ii = np.where(~dmag.mask)
			dmag_all[ref].append(np.array(bokmag-refmag)[ii])
			refclr_all[ref].append(np.array(refclr[ref][ii]))
		# iteratively fit a polynomial of increasing order to the
		# magnitude differences
		cterms = {}
		for ref in refs:
			_dmag = np.concatenate(dmag_all[ref]).flatten()
			_refclr = np.concatenate(refclr_all[ref]).flatten()
			mask = np.abs(_dmag) > 0.25
			for iternum in range(3):
				order = iternum+1
				tmp_dmag = np.ma.array(_dmag,mask=mask)
				fit = np.ma.polyfit(_refclr,tmp_dmag,order)
				magoff = sigma_clip(_dmag-np.polyval(fit,_refclr))
				mask = magoff.mask
			cterms[ref] = fit
			if savefit:
				_fit = fit.copy()
				# remove the lowest order term (absorbed by zeropoint)
				_fit[-1] = 0
				np.savetxt('config/bok2%s_%s_gicoeff.dat'%(ref,b),_fit)
		if doplots:
			for ref in refs:
				_dmag = np.concatenate(dmag_all[ref])
				_refclr = np.concatenate(refclr_all[ref])
				fig = plt.figure(figsize=(10,6))
				fig.subplots_adjust(0.09,0.1,0.98,0.92)
				ax1 = fig.add_subplot(211)
				#plt.scatter(_refclr,_dmag,
				#            s=3,c='0.2',edgecolor='none')
				ax1.hexbin(_refclr,_dmag,
				           cmap=plt.get_cmap('gray_r'),bins='log')
				ax1.axhline(0,c='c')
				xx = np.linspace(-1,5,100)
				ax1.plot(xx,np.polyval(cterms[ref],xx),c='r')
				try:
					oldcterms = np.loadtxt('config/bok2%s_%s_gicoeff.dat' %
					                       (ref,b))
					yoff = cterms[ref][-1]
					ax1.plot(xx,yoff+np.polyval(oldcterms,xx),c='orange')
				except:
					pass
				ax1.set_ylabel('%s(Bok) - %s(%s)'%(b,b,ref.upper()))
				ax1.set_ylim(-0.25,0.25)
				order = len(cterms[ref])-1
				clrstr = {'sdss':'gi','ps1':'gi'}[ref]
				polystr = ' '.join(['%+.5f*%s^%d'%(c,clrstr,order-d) 
				                      for d,c in enumerate(cterms[ref])])
				ax1.set_title(('$%s(Bok) - %s(%s) = '%(b,b,ref.upper())) +
				              polystr+'$',size=14)
				ax2 = fig.add_subplot(212)
				ax2.hexbin(_refclr,_dmag-np.polyval(cterms[ref],_refclr),
				           cmap=plt.get_cmap('gray_r'),bins='log')
				ax2.axhline(0,c='r')
				ax2.set_ylim(-0.15,0.15)
				for ax in [ax1,ax2]:
					if ref=='sdss':
						ax.axvline(0.4,c='b')
						ax.axvline(3.7,c='b')
						ax.set_xlabel('SDSS g-i')
						ax.set_xlim(-0.05,3.9)
					else:
						ax.axvline(0.4,c='b')
						ax.axvline(2.7,c='b')
						ax.set_xlabel('PS1 g-i')
						ax.set_xlim(-0.05,2.8)
					ax.axvline(0,c='MidnightBlue')
					ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
					ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
				fig.savefig('%s_%s_to_bok_gicolors.png'%(b,ref))

def calc_zeropoints(ref='ps1'):
	t = load_nov2015_data()
	refclr = {'sdss':t['sdss_mag_g']-t['sdss_mag_z'],
	           'ps1':t['ps1_mag_g']-t['ps1_mag_i']}[ref]
	zp = np.zeros((12,2,16))
	for j,b in enumerate('gr'):
		coeff = np.loadtxt('config/bok2%s_%s_coeff.dat'%(ref,b))
		refmag = t[ref+'_mag_'+b] + np.polyval(coeff,refclr.filled(0))
		refmag = refmag[:,np.newaxis]
		# find the zeropoint offsets for each image and amplifier independently
		for k,amp in enumerate(range(1,17)):
			is_amp = t['bok_amp_'+b]==amp
			for field in range(2):
				# split the two fields
				if field==0:
					is_field = t['sdss_ra'] > 345
				else:
					is_field = t['sdss_ra'] < 345
				mag = t['bok_mag_'+b].copy()
				mag.mask |= refclr.mask[:,np.newaxis]
				mag.mask[~(is_field[:,np.newaxis]&is_amp)] = True
				dmag = sigma_clip(refmag-mag,axis=0)
				if field==0:
					zp[:6,j,k] = dmag.mean(axis=0)
				else:
					zp[6:,j,k] = dmag.mean(axis=0)
	# np.savez('foo',zpEl=zp,zpADU=...,gain=)
	return zp

# TODO convert zps back to ADU and apply nominal ext corr



def dump_zeropoints(fieldname,band,byamp=False,showADU=True,**kwargs):
	from boketc import k_ext
	ampn = np.array(ampNums)
	tab = Table.read(fieldname+'_merged.fits')
	zpim,zpccd,zpamp = calc_zeropoints(tab,band,**kwargs)
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',
	                      'bokpipe_v0.2','nov15data')
	if showADU:
		obsdb = load_obsdb()
		fname = fieldname.replace('r_','bokr_')
		airmass = np.array([ obs['airmass'] for obs in obsdb
		                         if obs['objName'].startswith(fname) and
		                            obs['good'] ])
		zp0adu = zpim - 2.5*np.log10(1.375) - k_ext[band]*(airmass-1)
	for j in range(len(zpim)):
		print 'image %2d: ' % (j+1),
		for ccdNum in range(1,5):
			print '%6.3f ' % zpccd[j,ccdNum-1],
		if byamp:
			print
			for ccdNum in range(1,5):
				print '   ',
				for ai in range(4):
					print 'IM%-2d  %6.3f   ' % \
					   (ampn[ccdNum-1][ai],zpamp[j,ccdNum-1,ai]),
				print
		print '   %6.3f   %6.3f' % (zpim[j],zp0adu[j])
	print
	print ' average zeropoints: %6.3f   %6.3f' % (zpim.mean(),zp0adu.mean())
	if byamp and showADU:
		from bokpipe.bokproc import nominal_gain,ampOrder
		# to derive per-amp zeropoints using applied gain corrections
		gaindat = np.load(os.path.join(bokdir,'..','diagnostics',
		              'gainbal_20151112_%s.npz'%{'g':'g','r':'bokr'}[band]))
		gain = nominal_gain * np.product(gaindat['gainCor'],axis=0)
		# reorder arrays to match ampNumber
		ampgain = nominal_gain[np.array(ampOrder)-1][ampn-1]
		deltazp = zpim.mean() - zpamp
		ampdelta = deltazp.mean(axis=0)
		extcorr = (k_ext[band]*(airmass-1))[:,np.newaxis,np.newaxis]
		zp0adu_amp = zpamp - 2.5*np.log10(ampgain) - extcorr + ampdelta
		zp0adu_amp = zp0adu_amp.mean(axis=0)
		print
		print ' per-amplifier zeropoints in ADU:'
		for ccdNum in range(1,5):
			for ai in range(4):
				print '  IM%-2d  %6.3f  ' % \
				   (ampn[ccdNum-1][ai],zp0adu_amp[ccdNum-1,ai]),
			print
		print
		_zp = zp0adu_amp.flatten()
		_ii = _zp.argsort()
		print 'sorted by zeropoint(ADU):'
		for _i in _ii:
			print '  IM%-2d  %6.3f' % (ampn.flatten()[_i],_zp[_i])
		print ' median is %6.3f' % np.median(zp0adu_amp)
		print

def delta_zps_byamp(zpFit):
	for ccdi in range(4):
		ccdmean = np.mean(zpFit[:,ccdi,:])
		for ampj in range(4):
			ampmean = np.mean(zpFit[:,ccdi,ampj])
			print 'IM%d %8.4f' % (ampNums[ccdi][ampj],ccdmean-ampmean)

def dump_all_zeropoints():
	for fpfx,fsfx in [('s82cal','ra334'),('deep2','ra352')]:
		for filt in ['g','r']:
			fieldcatf = fpfx+filt+'_'+fsfx+'_merged.fits'
			fieldcat = Table.read(fieldcatf)
			print fieldcatf
			dump_zeropoints(fieldcat,filt[-1])
			print
		print

def check_scatter(_tab,band,j=None,colorcorr=False):
	from scipy.stats import scoreatpercentile
	piles = [25,50,75,90,95]
	for i in range(2):
		if i==0:
			tab = _tab[_tab['sdss_ra']<340]
			print '\nStripe 82 RA=334'
		else:
			tab = _tab[_tab['sdss_ra']>340]
			print '\nDeep2F3'
		print '%8s ' % 'mag',
		print ' -%2d%%- '*len(piles) % tuple(piles)
		for mag1 in np.arange(17,21.1,0.5):
			is_mag = ( (tab['sdss_mag_'+band]>(mag1-0.25)) & 
			           (tab['sdss_mag_'+band]<(mag1+0.25)) )
			if j is None:
				mag = tab['bok_meanmagcorr_'+band]
			else:
				mag = tab['bok_magcorr_'+band][:,j]
			ii = np.where(~mag.mask & is_mag)[0]
			dm = mag[ii] - tab['sdss_mag_'+band][ii]
			dist = scoreatpercentile(np.abs(dm-np.median(dm)),piles)
			print '%8.1f ' % mag1,
			print '%6d '*len(piles) % tuple(1e3*dist)
		print

def get_zeropoints(tab,zps,zptype):
	zpIm,zpCCD,zpAmp = zps
	if zptype=='image':
		zp = zpIm[:,np.newaxis]
	elif zptype=='ccd':
		# add zp=0 for ccd=0 (empty entries)
		zp = np.hstack([np.zeros((zpCCD.shape[0],1)),zpCCD])
		zp = np.choose(tab['ccdNum'],zp.transpose())[:,:,np.newaxis]
	elif zptype=='amp':
		ampIndex = get_amp_index(tab['X_IMAGE'],tab['Y_IMAGE'])
		ai = np.clip(4*(tab['ccdNum']-1) + ampIndex + 1, 0, 16)
		zp = np.hstack([np.zeros((zpAmp.shape[0],1)),zpAmp.reshape(-1,16)])
		zp = np.choose(ai,zp.transpose())[:,:,np.newaxis]
	return zp

def focalplanevar(tab,band,nbin=4,doplot=False,vr=0.015,shownum=False,
                  zpcorr=None):
	apNum = 2
	if zpcorr is None:
		zp = None
	else:
		zps = calc_zeropoints(tab,band,apNum=apNum)
		zp = get_zeropoints(tab,zps,zpcorr)
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
	s82r = Table.read('s82calr_ra334_merged.fits')
	plt.ioff()
	for filt,tab in [('g',s82g),('r',s82r)]:
		for zpcorr in ['image','ccd','amp']:
			fpim = focalplanevar(tab,filt,doplot=True,nbin=8,vr=0.015,
			                     zpcorr=zpcorr)
			plt.savefig('s82zp_%s_%s.png' % (filt,zpcorr))
	plt.ion()

def stripe82_seeing():
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'nov15data')
	for filt in 'gr':
		fields = ['s82cal%s_ra334_%s' % (filt,n) for n in '123456abcde']
		for field in fields:
			cat = fits.getdata(os.path.join(bokdir,filt,field+'.cat.fits'))
			ii = np.where((cat['FWHM_IMAGE']<10)&(cat['FLAGS']==0))[0]
			fwhm = np.ma.median(sigma_clip(cat['FWHM_IMAGE'][ii]))[0]
			print '%s %.2f %.2f' % (field,fwhm,fwhm*0.455)

def stripe82_aperphot_all(which='s82cal',redux='naoc',**kwargs):
	rdxdir = os.path.join(os.environ['BASSRDXDIR'],'reduced')
	if redux=='naoc':
		bokdir = os.path.join(rdxdir,'nov15data')
	elif redux=='naocv1':
		bokdir = os.path.join(rdxdir,'nov15data_OLD')
	elif redux=='idm':
		bokdir = os.path.join(rdxdir,'bokpipe_v0.2','nov15data')
	elif redux=='noao':
		bokdir = os.path.join(rdxdir,'noaocp','nov15data')
	else:
		raise ValueError
	aperRad = kwargs.pop('aperRad',7.0) / 0.455 # asec -> pix
	if which=='s82cal':
		s82cat = Table(load_stripe82_truth(ra_range=(332,336)),masked=True)
		nfields = 11
		fname = 's82cal%s_ra334_%s'
	else:
		s82cat = Table(load_stripe82_truth(ra_range=(349,354)),masked=True)
		nfields = 6
		fname = 'deep2%s_ra352_%s'
	exptime = np.array([100.]*6 + [25.,50.,100.,200.,400.])[:nfields]
	cols = ['x','y','counts','countsErr','peakCounts','flags','ccdNum']
	for k in cols:
		dtype = np.float32 if k not in ['flag','ccdNum'] else np.int32
		s82cat[k] = np.zeros((len(s82cat),nfields,2),dtype=dtype)
	for i in range(nfields):
		for j,filt in enumerate('gr'):
			field = fname % (filt,'123456abcde'[i])
			try:
				print 'calculating aperture phot for ',field,
				imf = os.path.join(bokdir,field+'.fits')
				bpmask = None #fitsio.FITS(imf.replace('.fits','.wht.fits'))
				ph = aper_phot_image(imf,s82cat['ra'],s82cat['dec'],[aperRad],
				                     badPixMask=bpmask,calc_peak=True)
				print
			except IOError:
				print '--> MISSING'
				continue
			for k in cols:
				s82cat[k][ph['objId'],i,j] = ph[k].squeeze()
	# objects with no coverage have ccdNum==0 in all frames/filters
	ii = np.where(s82cat['ccdNum'].sum(axis=-1).sum(axis=-1)>0)[0]
	s82cat = s82cat[ii]
	missing = s82cat['ccdNum']==0
	for k in cols:
		s82cat[k].mask |= missing
	s82cat['ampIndex'] = get_amp_index(s82cat['x'],s82cat['y'])
	s82cat['ampNum'] = 4*(s82cat['ccdNum']-1) + s82cat['ampIndex'] + 1
	if redux=='idm':
		# units of image are electrons
		s82cat['cps'] = s82cat['counts'] / exptime[np.newaxis,:,np.newaxis]
	else:
		# units of image are electrons/s
		s82cat['cps'] = s82cat['counts'].copy()
		s82cat['counts'] = s82cat['cps'] * exptime[np.newaxis,:,np.newaxis]
		s82cat['countsErr'] *= exptime[np.newaxis,:,np.newaxis]
	s82cat['meanCps'] = np.average(s82cat['cps'],weights=exptime,axis=1)
	s82cat['snr'] = s82cat['counts']/s82cat['countsErr']
	s82cat['nobs'] = (~s82cat['counts'].mask).sum(axis=1)
	ctsratio = np.ma.divide(s82cat['cps'],s82cat['meanCps'][:,np.newaxis,:])
	s82cat['dmag'] = -2.5*np.ma.log10(ctsratio)
	s82cat['sigMag'] = np.ma.std(s82cat['dmag'],axis=1)
	return s82cat

def selfcal(s82cat):
	exptime = np.array([100.]*6 + [25.,50.,100.,200.,400.])
	ii = np.where(np.any(s82cat['nobs']>6,axis=1))[0]
	imgcal = sigma_clip(s82cat['dmag'][ii],
	                    sigma=2.0,iters=3,axis=0).mean(axis=0)
	print imgcal
	cscl = 10**(-0.4*imgcal)
	s82cat['cpsCal'] = s82cat['cps']*cscl
	s82cat['meanCpsCal'] = np.average(s82cat['cpsCal'],
	                                  weights=exptime,axis=1)
	ctsratio = np.ma.divide(s82cat['cpsCal'],
	                        s82cat['meanCpsCal'][:,np.newaxis,:])
	s82cat['dmagCal'] = -2.5*np.ma.log10(ctsratio)
	s82cat['dmagCal'] -= sigma_clip(s82cat['dmagCal'][ii],
	             sigma=2.0,iters=3,axis=0).mean(axis=0).mean(axis=0)
	s82cat['sigMagCal'] = np.ma.std(s82cat['dmagCal'],axis=1)
	return s82cat

def sdsscal(s82cat,mode='amp'):
	exptime = np.array([100.]*6 + [25.,50.,100.,200.,400.])
	nimg = len(exptime)
	bokmag = -2.5*np.ma.log10(s82cat['cps']/exptime[:,np.newaxis])
	s82cat['dmagExt'] = 0*s82cat['dmag']
	s82cat['dmagExt'] = np.ma.masked
	gcterms = np.loadtxt('config/bok2sdss_g_gicoeff.dat')
	rcterms = np.loadtxt('config/bok2sdss_r_gicoeff.dat')
	gminusi = s82cat['psfmag_g']-s82cat['psfmag_i']
	sdssmag = {'g':s82cat['psfmag_g']+np.polyval(gcterms,gminusi),
	           'r':s82cat['psfmag_r']+np.polyval(rcterms,gminusi)}
	if mode=='amp':
		zp = np.ma.zeros((nimg,2,16))
		for amp in range(1,17):
			isamp = s82cat['ampNum'] == amp
			for j in range(nimg):
				for k in range(2):
					b = 'gr'[k]
					good = s82cat['nobs'][:,k] >= 1
					ii = np.where(good & isamp[:,j,k])[0]
					if len(ii)==0:
						zp[j,k,amp-1] = np.ma.masked
						continue
					dmag = np.ma.subtract(bokmag[ii,j,k],sdssmag[b][ii])
					zp[j,k,amp-1] = np.ma.median(sigma_clip(dmag,sigma=2.0))
					s82cat['dmagExt'][ii,j,k] = dmag - zp[j,k,amp-1]
					s82cat['dmagExt'][ii,j,k].mask |= dmag.mask
	elif mode=='ccd':
		zp = np.ma.zeros((nimg,2,4))
		for ccd in range(1,5):
			isccd = s82cat['ccdNum'] == ccd
			for j in range(nimg):
				for k in range(2):
					b = 'gr'[k]
					good = s82cat['nobs'][:,k] >= 1
					ii = np.where(good & isccd[:,j,k])[0]
					if len(ii)==0:
						zp[j,k,ccd-1] = np.ma.masked
						continue
					dmag = np.ma.subtract(bokmag[ii,j,k],sdssmag[b][ii])
					zp[j,k,ccd-1] = np.ma.median(sigma_clip(dmag,sigma=2.0))
					s82cat['dmagExt'][ii,j,k] = dmag - zp[j,k,ccd-1]
					s82cat['dmagExt'][ii,j,k].mask |= dmag.mask
	elif mode=='image':
		raise ValueError
	else:
		raise ValueError
	s82cat['sigMagExt'] = np.ma.std(s82cat['dmagExt'],axis=1)
	return s82cat,zp

def focalplanevar2(s82cat,band,nbin=4,frameno=None,dmk='dmag',
                   doplot=False,vr=0.015,shownum=False):
	tab = s82cat
	bk = 'gr'.find(band)
	nx = 4096 // nbin
	ny = 4032 // nbin
	if frameno is None:
		s = np.s_[:,:,bk]
	else:
		s = np.s_[:,frameno,bk]
	xi = (tab['x'][s]/nx).astype(np.int32)
	yi = (tab['y'][s]/ny).astype(np.int32)
	fpIm = np.zeros((4,nbin,nbin))
	for ccdi in range(4):
		for i in range(nbin):
			for j in range(nbin):
				ii = np.where((tab['ccdNum'][s]==ccdi+1)&(yi==i)&(xi==j))
				dm = sigma_clip(s82cat[dmk][s][ii])
				if shownum:
					fpIm[ccdi,i,j] = (~dm.mask).sum()
				else:
					fpIm[ccdi,i,j] = np.ma.mean(dm)
	if doplot:
		if vr is None:
			vmin,vmax = None,None
		else:
			try:
				vmin,vmax = vr
			except:
				vmin,vmax = -vr,vr
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
				                vmin=vmin,vmax=vmax,cmap=plt.cm.RdBu)
			plt.title('CCD%d'%(ccdi+1))
			if pnum==1:
				cbax = fig.add_axes([0.1,0.98,0.8,0.015])
				cb = fig.colorbar(_im,cax=cbax,orientation='horizontal')
				#if not shownum:
				#	cb.locator = ticker.MultipleLocator(0.005)
				#cb.update_ticks()
			ax.xaxis.set_ticks([])
			ax.yaxis.set_ticks([])
	return fpIm

def stripe82_linearity_plot(s82tab,filt='g',peak=False):
	from scipy.stats import scoreatpercentile
	from astrotools.idmstuff import binmean
	countsk = 'counts' if not peak else 'peakCounts'
	exptime = np.array([25.,50.,100.,200.,400.])
	magrange = (s82tab['psfmag_'+filt]>17) & (s82tab['psfmag_'+filt]<22)
	plt.figure(figsize=(8.5,9.5))
	plt.subplots_adjust(0.10,0.06,0.96,0.98,0.02,0.02)
	bk = 'gr'.find(filt)
	ii = np.where(s82tab[countsk].mask[:,-5:,bk].sum(axis=1)==0)[0]
	counts = s82tab[countsk][ii,-5:,bk].filled(-1)
	cps = s82tab['cps'][ii,-5:,bk].filled(-1)
	mean_cps = s82tab['meanCps'][ii,bk].filled(-1)
	amps = s82tab['ampNum'][ii,-1,bk].filled(-1) # just use last im
	flags = s82tab['flags'][ii,-5:,bk].filled(9999)
	magrange = magrange[ii]
	expected_cts = mean_cps[:,np.newaxis] * exptime
	for ampNum in range(1,17):
		ax = plt.subplot(8,2,ampNum)
		ii = np.where((amps==ampNum) & magrange & (mean_cps>0) &
		              np.all(flags==0,axis=1))[0]
		xall,yall = [],[]
		for j in range(len(exptime)):
			ctsratio = counts[ii,j]/expected_cts[ii,j]
			plt.scatter(np.log10(mean_cps[ii]),ctsratio,
			            s=7,c='gray',edgecolor='none')
			xall.append(np.log10(mean_cps[ii]))
			yall.append(ctsratio)
		xall = np.concatenate(xall)
		yall = np.concatenate(yall)
		plt.axhline(1.0,c='r')
		xr = scoreatpercentile(xall,[5,99])
		xb,yb,yv = binmean(xall,yall,
		                   np.linspace(xr[0],xr[1],10),std=True,median=True,
		                   clip=True)
		plt.errorbar(xb,yb,yv,fmt='ks',ms=4)
		if countsk=='counts':
			plt.xlim(1.4,3.4)
			plt.ylim(0.921,1.079)
		else:
			plt.xlim(1.0,2.4)
			plt.ylim(0.89,1.11)
		ax.yaxis.set_major_locator(ticker.MultipleLocator(0.04))
		ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
		if ampNum % 2 == 0:
			ax.yaxis.set_ticks([])
		if ampNum < 15:
			ax.xaxis.set_ticks([])
		plt.text(3.15,1.05,'IM%d'%ampNum)
	plt.figtext(0.5,0.01,'$log(<flux>)$',size=14,ha='center')
	plt.figtext(0.01,0.5,r'$flux / <flux>$',size=14,va='center',
	            rotation='vertical')

def rename_proc_files(which='naoc'):
	rdxdir = os.environ['BASSRDXDIR']
	if which=='naoc':
		bokdir = os.path.join(rdxdir,'reduced','20151111')
		outdir = os.path.join(rdxdir,'reduced','nov15data')
	elif which=='noao':
		bokdir = os.path.join(rdxdir,'BOK_CP','CP20151111V0')
		outdir = os.path.join(rdxdir,'reduced','noaocp','nov15data')
	else:
		raise ValueError
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	log = Table.read('bass_Nov2015toFeb2016.fits')
	ii = np.where(log['utDate']=='20151112')[0]
	for i in ii:
		for sfx in ['','.wht']:
			field = log['objName'][i]
			if not (field.startswith('cosmos') or field.startswith('deep2')
			         or field.startswith('s82cal')):
				continue
			outfn = os.path.join(outdir,field.replace('bokr','r')+sfx+'.fits')
			if os.path.exists(outfn): continue
			filt = log['filter'][i]
			if which=='naoc' and True:
				filt = filt.replace('bokr','r')
			if which=='noao':
				nsfx = 'oow' if sfx=='.wht' else 'ooi'
				_fn = log['fileName'][i].replace('ori',nsfx)+'_%s_v0' % filt
				infn = os.path.join(bokdir,_fn)+'.fits.fz'
				os.system('funpack -O %s %s' % (outfn,infn))
				continue
			# for NAOC images combine 4 CCDs into single MEF file
			d = log['DTACQNAM'][i].split('.')
			hdus = [fits.PrimaryHDU()]
			try:
				for ccd in range(1,5):
					fn = 'p'+d[0][1:]+filt+d[1]+'_%d%s.fits' % (ccd,sfx)
					im,hdr = fits.getdata(os.path.join(bokdir,fn),header=True)
					for k in ['CTYPE1','CTYPE2']:
						hdr[k] = hdr[k].replace('TAN','TPV')
					hdus.append(fits.ImageHDU(im,hdr,'CCD%d'%ccd))
				hdul = fits.HDUList(hdus)
				hdul.writeto(outfn)
			except IOError:
				pass

def plot_pointings():
	from matplotlib.patches import Rectangle
	from astrotools.idmstuff import radec_fromstr
	filt = 'g'
	#fields = ['cosmos%s_ra150_%s' % (filt,n) for n in '123456']
	fields = ['deep2%s_ra352_%s' % (filt,n) for n in '123456']
	_minmax = lambda a,b: (a,b) if a<=b else (b,a)
	plt.figure()
	ax = plt.subplot(111)
	for f in fields:
		imf = os.path.join(os.environ['BASSRDXDIR'],'reduced','bokpipe_v0.2',
		                   'nov15data',filt,f+'.fits')
		hdus = fits.open(imf)
		for hdu in hdus[1:]:
			w = WCS(hdu.header)
			ra0,dec0 = w.all_pix2world(1,1,1,ra_dec_order=True)
			ra1,dec1 = w.all_pix2world(4096,4032,1,ra_dec_order=True)
			ra0,ra1 = _minmax(ra0,ra1)
			dec0,dec1 = _minmax(dec0,dec1)
			rect = Rectangle((ra0,dec0),ra1-ra0,dec1-dec0,
			                 alpha=0.5,color='0.2')
			ax.add_patch(rect)
	moscoo = ['10:00:28.00  02:12:20.99', '10:00:19.06  02:17:11.00',
	          '10:00:25.80  02:15:56.99', '10:00:21.25  02:13:35.00',
	          '10:00:23.52  02:14:45.99']
	for c in moscoo:
		ra,dec = radec_fromstr(c,True)
		#plt.scatter(ra,dec,c='r',s=100,lw=2,marker='x')
		rect = Rectangle((ra-0.3,dec-0.3),0.6,0.6,
		                 color='r',fill=False)
		ax.add_patch(rect)
	d2f3coo = ['23:25:00 -00:05:00', '23:35:00 -00:05:00',
	           '23:35:00 +00:27:00', '23:25:00 +00:27:00',
	           '23:25:00 -00:05:00']
	d2f3 = np.array([radec_fromstr(c,True) for c in d2f3coo])
	print d2f3.shape
	plt.plot(d2f3[:,0],d2f3[:,1],c='r')
	#plt.xlim(149.49,151.3)
	#plt.ylim(1.05,3.05)
	plt.xlim(351,353.5)
	plt.ylim(-1.4,1.0)

def get_colors():
	apNum = 2
	pfx = 's82cal%s_ra334'
	mag = {}
	for b in 'gr':
		tab = Table.read(pfx%b+'_merged.fits')
		zps = calc_zeropoints(tab,b,apNum=apNum)
		zp = get_zeropoints(tab,zps,'amp')
		aperMag = flux2mag(tab,b,'APER',zp=zp)[:,:,apNum]
		meanMag = sigma_clip(aperMag,axis=1).mean(axis=1)
		refMag = tab['tMag']
		ii = np.where(( tab['tMag']>16.5) & (tab['tMag']<21.5) &
		                ~meanMag.mask )[0]
		mag[b] = {'bok_'+b:meanMag[ii].filled(0),'ref_'+b:refMag[ii],
		          'idx':tab['tIndex'][ii]}
	tab = join(Table(mag['g']),Table(mag['r']),keys='idx')
	tab['bok_gr'] = tab['bok_g'] - tab['bok_r']
	tab['ref_gr'] = tab['ref_g'] - tab['ref_r']
	return tab

def color_terms(tab=None,nIter=3):
	if not tab:
		tab = get_colors()
	plt.figure(figsize=(9,4))
	plt.subplots_adjust(0.08,0.12,0.95,0.95,0.3,0.3)
	xedges = np.arange(-0.5,2.01,0.05)
	yedges = np.arange(-0.5,0.51,0.02)
	xbins = xedges[:-1] + np.diff(xedges)/2
	ybins = yedges[:-1] + np.diff(yedges)/2
	cfits = []
	for pnum,b in enumerate('gr',start=1):
		dmag = tab['bok_'+b] - tab['ref_'+b]
		dmag = sigma_clip(dmag,sigma=5.0)
		ii = np.where(tab['ref_gr']>1.3)[0]
		dmag[ii] = np.ma.masked
		print len(dmag),dmag.mask.sum()
		for iterNum in range(nIter):
			fit = np.ma.polyfit(tab['ref_gr'],dmag,1)
			resid = dmag - np.polyval(fit,tab['ref_gr'])
			dmag[np.abs(resid)>3.5*resid.std()] = np.ma.masked
			print iterNum,fit,dmag.mask.sum()
		print
		cfits.append(fit)
		xx = np.array([-0.5,2.0])
		n,_,_ = np.histogram2d(tab['ref_gr'],dmag,[xedges,yedges])
		ax = plt.subplot(1,2,pnum)
		plt.axhline(0,c='m',ls='--')
		plt.scatter(tab['ref_gr'],dmag.data,
		            edgecolor='none',alpha=0.7,s=1,c='k')
		plt.contour(xbins,ybins,n.transpose(),colors='b')
		plt.plot(xx,np.polyval(fit,xx),c='r',lw=2)
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
		ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
		plt.xlim(0.1,1.7)
		plt.ylim(-0.35,0.25)
		plt.ylabel('Bok %s - SDSS %s'%(b,b),size=11)
		plt.xlabel('SDSS g-r',size=11)
		plt.text(0.2,-0.3,
		         r'$%s(Bok) = %s(SDSS) %.2f\times(g-r) + %.2f$' %
		         ((b,b)+tuple(fit)),
		         size=13,bbox=dict(facecolor='w',alpha=0.8))
	np.savetxt('stripe82cterms.dat',cfits)


def etc_check():
	from boketc import texp_onsky,k_ext,bok_zpt0
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'nov15data')
	ebv = 0.04 # rough average in the field
	exptime = 100. # all the same
#	deep2g = Table.read('deep2g_ra352_merged.fits')
#	deep2r = Table.read('deep2r_ra352_merged.fits')
#	zps = {b:get_zeropoints(tab,calc_zeropoints(tab,b),'image').squeeze()
#	          for b,tab in [('g',deep2g),('r',deep2r)]}
	zps = {'g':[25.909,25.914,25.922,25.913,25.921,25.926,25.934,
	            25.932,25.936,25.936,25.938],
	       'r':[25.753,25.761,25.756,25.750,25.707,25.732,25.739,
	            25.746,25.749,25.737,25.753]}
	for filt in 'gr':
		print '[%s] %10s  %4s %6s %6s %5s %7s %10s' % \
		        (filt,'field  ','secz','skyext','skymag','fwhm','t_est',
		         'depth_fac')
		fields = ['deep2%s_ra352_%s' % (filt,n) for n in '123456']
		for field,zp in zip(fields,zps[filt]):
			imf = os.path.join(bokdir,filt,field+'.fits')
			hdr0 = fits.getheader(imf,0)
			hdr1 = fits.getheader(imf,1)
			airmass = hdr0['airmass']
			b = {'g':'g','r':'bokr'}[filt]
			zp += -2.5*np.log10(1.375) - k_ext[b]*(airmass-1)
			skyextinction = np.clip(bok_zpt0[b]-zp,0,np.inf)
			cat = fits.getdata(os.path.join(bokdir,filt,field+'.cat.fits'))
			ii = np.where((cat['MAG_AUTO']-cat['MAG_PSF']>-0.08) &
			              (cat['MAGERR_AUTO']<0.03) &
			              (cat['FLAGS']==0))[0]
			fwhm = np.ma.median(sigma_clip(cat['FWHM_IMAGE'][ii]))[0]
			fwhm *= 0.455
			skyeps = hdr1['SKYVAL']/exptime
			skymag = -2.5*np.log10(skyeps) + bok_zpt0[b]
			t = texp_onsky(b,airmass,ebv,skyextinction,
			               skyeps=skyeps,fwhm=fwhm)
			dfac = -2.5*0.5*np.log10(100/(t/3))
			print '%s  %4.2f %6.2f %6.2f %5.2f %7.1f %10.2f' % \
			        (field,airmass,skyextinction,skymag,fwhm,t/3,dfac)
		print

def compare_scatter(phot1,phot2,names,minnobs=4,sigk='sigMag'):
	m1,m2 = srcor(phot1['ra'],phot1['dec'],phot2['ra'],phot2['dec'],1.0)
	plt.figure(figsize=(10,5))
	plt.subplots_adjust(0.08,0.12,0.95,0.92,0.25)
	for j,b in enumerate('gr'):
		ii = np.where((phot1['nobs'][m1,j]>=minnobs) &
		              (phot2['nobs'][m2,j]>=minnobs))[0]
		plt.subplot(1,2,j+1)
		for p,c,l in zip([phot1[m1[ii]],phot2[m2[ii]]],'bg',names):
			ii = np.where(p[sigk][:,j] > 0)[0]
			plt.scatter(p['psfmag_'+b][ii],p[sigk][ii,j],s=7,c=c,
			            edgecolor='none',alpha=0.7,label=l)
		plt.title('%s band' % b)
		plt.legend(loc='upper left')
		plt.xlim(15.7,{'g':22.5,'r':22.1}[b])
		plt.ylim(-0.005,0.35)
		plt.xlabel('SDSS %s mag' % b)
		plt.ylabel(r'$\sigma(%s)$' % b)

