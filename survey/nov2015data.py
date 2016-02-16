#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table,join
from astropy.stats import sigma_clip
from astropy.wcs import WCS

import matplotlib.pyplot as plt
from matplotlib import ticker

from boketc import bok_zpt0

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

def match_stripe82():
	from collections import defaultdict
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',
	                      'bokpipe_v0.2','nov15data')
	s82 = load_stripe82_truth()
	cols = ['NUMBER','FLUX_APER','FLUXERR_APER','FLUX_AUTO','FLUXERR_AUTO',
	        'BACKGROUND','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000',
	        'ELONGATION','ELLIPTICITY','FWHM_IMAGE','FLAGS','CLASS_STAR',
	        'FLUX_RADIUS','FLUX_MAX','SNR_WIN','FLUX_PSF','FLUXERR_PSF']
	nAper = 3
	cdtypes = {'NUMBER':'i4','ALPHA_J2000':'f8','DELTA_J2000':'f8',
	           'FLAGS':'i8'}
	for filt in ['g','r']:
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
			fn = fields[0][:-2].replace('bokr','r')
			tab.write(fn+'_merged.fits',overwrite=True)

def flux2mag(tab,band,which='PSF',zp=None):
	assert which in ['APER','AUTO','PSF']
	if zp is None:
		zp = bok_zpt0[{'g':'g','r':'bokr'}[band]]
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

def get_amp_index(x,y):
	nx = 4096 // 2
	ny = 4032 // 2
	xi = (x/nx).astype(np.int32)
	yi = (y/ny).astype(np.int32)
	ampIndex = 2*yi + xi
	return ampIndex

def calc_zeropoints(tab,band,apNum=2,zp=None,savedelta=False):
	if zp is None:
		zp = bok_zpt0[{'g':'g','r':'bokr'}[band]]
	aperMag = flux2mag(tab,band,'APER',zp=zp)
	ref_star = (tab['tMag']>16.5) & (tab['tMag']<19.7)
	nImages = tab['NUMBER'].shape[1]
	ampIndex = get_amp_index(tab['X_IMAGE'],tab['Y_IMAGE'])
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

def dump_zeropoints(fieldname,band,byamp=False,showADU=True,**kwargs):
	from boketc import k_ext
	#from bokpipe.bokproc import nominal_gain
	tab = Table.read(fieldname+'_merged.fits')
	zpim,zpccd,zpamp = calc_zeropoints(tab,band,**kwargs)
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',
	                      'bokpipe_v0.2','nov15data')
	# to derive per-amp zeropoints using applied gain corrections (not used)
	#gaindat = np.load(os.path.join(bokdir,'..','diagnostics',
	#              'gainbal_20151112_%s.npz'%{'g':'g','r':'bokr'}[band]))
	#gain = nominal_gain * np.product(gaindat['gainCor'],axis=0)
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
			if not byamp:
				continue
			for ai in range(4):
				print '  IM%d  %6.3f  ' % \
				   (ampNums[ccdNum-1][ai],zpamp[j,ccdNum-1,ai]),
			print
		print '   %6.3f   %6.3f' % (zpim[j],zp0adu[j])
	print
	print ' average zeropoints: %6.3f   %6.3f' % (zpim.mean(),zp0adu.mean())

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

def stripe82_phot(imageFile,s82cat,aperRad=2.5):
	from bokpipe.bokphot import aper_phot_image
	aperRad /= 0.455
	ph = aper_phot_image(imageFile,s82cat['ra'],s82cat['dec'],[aperRad])
	ph['refMag'] = s82cat['psfMag_g'][ph['objId']]
	ph['ra'] = s82cat['ra'][ph['objId']]
	ph['dec'] = s82cat['dec'][ph['objId']]
	return ph

def stripe82_linearity(filt,**kwargs):
	bokdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',
	                      'bokpipe_v0.2','nov15data')
	s82all = load_stripe82_truth(ra_range=(332,336))
	fields = ['s82cal%s_ra334_%s' % (filt,n) for n in 'abcde']
	ph = [ stripe82_phot(os.path.join(bokdir,filt,field+'.fits'),
	                     s82all,**kwargs)
	                 for field in fields ]
	tab = np.zeros(len(s82all),
	               dtype=[('x','5f4'),('y','5f4'),
	                      ('counts','5f4'),('countsErr','5f4'),
	                      ('flags','5i4'),('ccdNum','i4')])
	for j,p in enumerate(ph):
		for k in ['x','y','counts','countsErr','flags']:
			tab[k][p['objId'],j] = p[k].squeeze()
		tab['ccdNum'][p['objId']] = p['ccdNum']
	# actually had coverage
	ii = np.where(tab['ccdNum']>0)[0]
	tab = tab[ii]
	tab = Table(tab)
	tab['refMag'] = s82all['psfMag_'+filt[-1]][ii]
	tab['refErr'] = s82all['psfMagErr_'+filt[-1]][ii]
	tab['ampIndex'] = get_amp_index(tab['x'][:,-1],tab['y'][:,-1])
	tab['ampNum'] = 4*(tab['ccdNum']-1) + tab['ampIndex'] + 1
	return tab

def stripe82_linearity_plot(s82tab):
	m = 19.5
	dm = 2.5
	exptime = np.array([25.,50.,100.,200.,400.])
	plt.figure(figsize=(7.5,9.5))
	if True:
		magrange = (s82tab['refMag']>m-dm) & (s82tab['refMag']<m+dm)
		for ampNum in range(1,17):
			plt.subplot(8,2,ampNum)
			ii = np.where((s82tab['ampNum']==ampNum) & magrange &
#			              (s82tab['counts'][:,-1]>0) &
			              np.all(s82tab['flags']==0,axis=1))[0]
			cps = s82tab['counts'][ii] / exptime
			mean_cps = np.average(cps,weights=exptime,axis=-1)
			for j in range(5):
				#xx = exptime[j] + 10*(np.random.rand(len(ii))-0.5)
				xx = mean_cps
				expected_cts = mean_cps * exptime[j]
				#ctsratio = s82tab['counts'][ii,j]/s82tab['counts'][ii,-1]
				#ctsratio *= exptime[-1]/exptime[j]
				ctsratio = s82tab['counts'][ii,j]/expected_cts
				plt.scatter(xx,ctsratio,
				            s=10,c='gray',edgecolor='none')
#				plt.errorbar(exptime[j],ctsratio.mean(),ctsratio.std(),
#				             fmt='bs')
			plt.axhline(1.0,c='r')
			plt.ylim(0.96,1.04)
			plt.xscale('log')
			plt.xlim(20,8e3)

def plot_pointings():
	from matplotlib.patches import Rectangle
	from astrotools.idmstuff import radec_fromstr
	filt = 'g'
	fields = ['cosmos%s_ra150_%s' % (filt,n) for n in '123456']
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
	plt.xlim(149.49,151.3)
	plt.ylim(1.05,3.05)

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

def color_terms(tab,nIter=3):
	plt.figure(figsize=(9,4))
	plt.subplots_adjust(0.08,0.12,0.95,0.95,0.3,0.3)
	xedges = np.arange(-0.5,2.01,0.05)
	yedges = np.arange(-0.5,0.51,0.02)
	xbins = xedges[:-1] + np.diff(xedges)/2
	ybins = yedges[:-1] + np.diff(yedges)/2
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

