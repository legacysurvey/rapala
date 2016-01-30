#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import astropy.io.ascii as ascii_io
import fitsio

import bass
import bokextract

datadir = '/global/scratch2/sd/imcgreer/'
ndwfs_starfile = datadir+'ndwfs/starcat.fits'
bootes_sdss_starfile = datadir+'ndwfs/sdss_bootes_gstars.fits'
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

def panstarrs_md_tiles(observed=True):
	tiles = {}
	for field,ra,dec in [('MD03',130.592,+44.317),
                         ('MD05',161.917,+58.083),
                         ('MD06',185.000,+47.117),
                         ('MD07',213.704,+53.083),
                         ('MD08',242.787,+54.950)]:
		dra = 3.5/np.cos(np.radians(dec))
		tiles[field] = bass.region_tiles(ra-dra,ra+dra,dec-3.5,dec+3.5,
		                                 observed=observed)
	return tiles

def check_fields_list():
	files = [ t['utDate']+'/'+t['fileName']+'.fits.gz'
	                 for tiles in [cfhtw3_tiles(),ndwfs_tiles()] 
	                      for t in tiles ]
	with open('checkfields_tiles.txt','w') as f:
		f.write('\n'.join(sorted(files)))

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

def srcorXY(x1,y1,x2,y2,maxrad):
	sep = np.sqrt( (x1[:,np.newaxis]-x2[np.newaxis,:])**2 + 
	               (y1[:,np.newaxis]-y2[np.newaxis,:])**2 )
	ii = sep.argmin(axis=1)
	m1 = np.arange(len(x1))
	jj = np.where(sep[m1,ii] < maxrad)[0]
	return m1[jj],ii[jj]

def match_objects(objs,tiles):
	objpars = [('g_number','f4'),('g_ra','f8'),('g_dec','f8'),
	           ('g_x','f4'),('g_y','f4'),
	           ('g_autoMag','f4'),('g_autoMagErr','f4'),
	           ('g_autoFlux','f4'),('g_autoFluxErr','f4'),
	           ('g_psfMag','f4'),('g_psfMagErr','f4'),
	           ('g_psfFlux','f4'),('g_psfFluxErr','f4'),
	           ('g_elongation','f4'),('g_ellipticity','f4'),
	           ('g_flags','i4'),('g_fluxRad','f4')]
	tilepars = [('g_utDate','S8'),('g_expTime','f4'),
	            ('g_tileId','i4'),('g_ditherId','i4'),('g_ccdNum','i4')]
	dtype = objs.dtype.descr + objpars + tilepars
	skeys = ['NUMBER','ALPHA_J2000','DELTA_J2000','X_IMAGE','Y_IMAGE',
	         'MAG_AUTO','MAGERR_AUTO','FLUX_AUTO','FLUXERR_AUTO',
	         'MAG_PSF','MAGERR_PSF','FLUX_PSF','FLUXERR_PSF',
	         'ELONGATION','ELLIPTICITY',
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
			ii = np.where( (objs['ra']>cat['ALPHA_J2000'].min()+3e-3) &
			               (objs['ra']<cat['ALPHA_J2000'].max()-3e-3) &
			               (objs['dec']>cat['DELTA_J2000'].min()+3e-3) &
			               (objs['dec']<cat['DELTA_J2000'].max()-3e-3) )[0]
			if len(ii)==0:
				continue
			m1,m2 = srcor(objs['ra'][ii],objs['dec'][ii],
			              cat['ALPHA_J2000'],cat['DELTA_J2000'],2.5)
			print '  ccd%d %d/%d' % (ccdNum,len(m1),len(ii)),
			matches.extend( [ tuple(objs[i]) +
			                  tuple([cat[k][j] for k in skeys]) +
			                  tuple([t[k] for k in tkeys]) + (ccdNum,)
			                     for i,j in zip(ii[m1],m2) ] )
			uu = np.delete(np.arange(len(ii)),m1)
			matches.extend( [ tuple(objs[i]) +
			                  tuple([0]*len(skeys)) + 
			                  tuple([t[k] for k in tkeys]) + (ccdNum,)
			                     for i in ii[uu] ] )
		print
	matches = np.array(matches,dtype=dtype)
	print 'finished with ',matches.size
	return matches

def depth_plots(matches,g_ref,gname,bypriority=True,aper='psf',**kwargs):
	assert aper in ['psf','auto']
	fluxk = 'g_%sFlux' % aper
	errk = 'g_%sFluxErr' % aper
	#
	m = np.where( (matches[fluxk]>0) & (matches[errk]>0) )[0]
	gSNR = matches[fluxk][m] / matches[errk][m]
	if bypriority:
		fig1 = plt.figure(figsize=(10,8))
		plt.subplots_adjust(0.07,0.07,0.97,0.96,0.27,0.27)
	else:
		fig1 = plt.figure(figsize=(5,4.5))
		plt.subplots_adjust(0.13,0.12,0.98,0.94)
	for i in range(4):
		if bypriority:
			ax = plt.subplot(2,2,i+1)
		else:
			if i>0: break
			ax = plt.subplot(1,1,i+1)
		if i==0:
			ii = np.where(matches['g_ditherId'][m] > 0)[0]
		else:
			ii = np.where(matches['g_ditherId'][m] == i)[0]
		ax.hexbin(g_ref[m[ii]],np.log10(gSNR[ii]),
		          bins='log',cmap=plt.cm.Blues)
		ax.axhline(np.log10(5.0),c='r',lw=1.3,alpha=0.7)
		ax.plot([24.0-2.5*np.log10(np.sqrt(3))]*2,np.log10([3,8]),c='m',lw=1.5)
		ax.set_xlim(17.2,24.5)
		ax.set_ylim(np.log10(2),np.log10(500))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
		ax.yaxis.set_major_locator(ticker.FixedLocator(np.log10(
		      [2,5,10,20,50,100,200])))
		ax.yaxis.set_major_formatter(ticker.FuncFormatter(
		      lambda x,pos: '%d' % np.round(10**x)))
		ax.set_xlabel(gname+' mag')
		ax.set_ylabel('BASS %s flux/err' % aper.upper())
		if i==0:
			ax.set_title('all tiles')
		else:
			ax.set_title('P%d tiles' % i)
	#
	mbins = np.arange(18.,24.01,0.1)
	fig2 = plt.figure(figsize=(8,4))
	plt.subplots_adjust(0.07,0.14,0.97,0.97,0.25)
	ax1 = plt.subplot(121)
	ax2 = plt.subplot(122)
	for i in range(4):
		if i==0:
			ii = np.where(matches['g_ditherId'] > 0)[0]
		else:
			if not bypriority: break
			ii = np.where(matches['g_ditherId'] == i)[0]
		jj = np.where(matches[errk][ii]>0)[0]
		g5sig = ( matches[fluxk][ii[jj]] / matches[errk][ii[jj]] ) > 5.0
		tot,_ = np.histogram(g_ref[ii],mbins)
		det,_ = np.histogram(g_ref[ii[jj]],mbins)
		det5,_ = np.histogram(g_ref[ii[jj[g5sig]]],mbins)
		ax1.plot(mbins[:-1],det.astype(np.float)/tot,drawstyle='steps-pre',
		         c=['black','blue','green','DarkCyan'][i],lw=1.3,
		         label=['all','P1','P2','P3'][i])
		ax2.plot(mbins[:-1],det5.astype(np.float)/tot,drawstyle='steps-pre',
		         c=['black','blue','green','DarkCyan'][i],lw=1.3,
		         label=['all','P1','P2','P3'][i])
	ax1.set_xlabel(gname+' mag')
	ax2.set_xlabel(gname+' mag')
	ax1.set_ylabel('fraction detected')
	ax2.set_ylabel('fraction detected 5 sig')
	ax1.legend(loc='lower left')
	if kwargs.get('saveplots',False):
		figname = kwargs.get('figname','blah')
		figext = kwargs.get('figtype','png')
		fig1.savefig(figname+'_depth.'+figext)
		fig2.savefig(figname+'_complete.'+figext)



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
		w = rfits[1].where('FWHM_IMAGE < 7 && MAG_AUTO < 24.0 && FLAGS == 0')
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

def ndwfs_depth(**kwargs):
	kwargs.setdefault('figname','ndwfs')
	ndwfsm = fitsio.read('ndwfs_match.fits')
	Bw = ndwfsm['autoMag'][:,0]
	Bw_minus_R = ndwfsm['autoMag'][:,0] - ndwfsm['autoMag'][:,1]
	NDWFSg = np.choose(Bw_minus_R <= 1.45, 
	                   [ Bw - (0.23*Bw_minus_R + 0.25),
	                     Bw - (0.38*Bw_minus_R + 0.05) ])
	#
	m = np.where( np.all(ndwfsm['autoMag'][:,:2]> 0,axis=1) &
	              np.all(ndwfsm['autoMag'][:,:2]<30,axis=1) )[0]
	depth_plots(ndwfsm[m],NDWFSg[m],'NDWFS g-ish',**kwargs)



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

def cfhtls_depth(**kwargs):
	kwargs.setdefault('figname','cfhtls')
	cfhtlsm = fitsio.read('cfhtlswide_match.fits')
	m = np.where( (cfhtlsm['psfMag'][:,1]> 0) &
	              (cfhtlsm['psfMag'][:,1]<30) )[0]
	depth_plots(cfhtlsm[m],cfhtlsm['psfMag'][m,1],'CFHTLS g',bypriority=False,
	            **kwargs)

bok_gain_2015 = [ 1.3325, 1.5225, 1.415, 1.47 ]
bok_rn_2015 = [ 7.94, 9.54, 11.81, 8.91 ]

def cfhtls_depth_compare():
	import itertools
	import boketc
	import bokdepth
	tiles = cfhtw3_tiles(observed=True)
	cfhtlsm = fitsio.read('stuff/cfhtlswide_match.fits')
	m = np.where( (cfhtlsm['psfMag'][:,1]>20) &
	              (cfhtlsm['psfMag'][:,1]<30) )[0]
	m = cfhtlsm[m]
	for ccdNum in range(1,5):
		ents = []
		for ti,t in enumerate(tiles):
			print ccdNum,ti,len(tiles)
			ii = np.where( (m['g_tileId'] == t['tileId']) &
			               (m['g_ditherId'] == t['ditherId']) &
			               (m['g_ccdNum'] == ccdNum) &
			               (m['g_psfFlux'] != 0) )[0]
			if len(ii)==0:
				continue
			impath = os.path.join(bass.rdxdir,t['utDate'],'ccdproc3',
			                      t['fileName']+'_ccd%d_pv.fits'%ccdNum)
			psfpath = os.path.join(bass.rdxdir,t['utDate'],'ccdproc3',
			                       t['fileName']+'_ccd%d.ldac_cat.psf'%ccdNum)
			if not os.path.exists(impath):
				print ' ... %s does not exist, skipping' % impath
				continue
			gain = bok_gain_2015[ccdNum-1]
			rdnoise = bok_rn_2015[ccdNum-1]
			rmsADU,rmsEl,A,skyADU = bokdepth.calc_processed_image_rms(
			                                       impath,psfpath,
			                                       gain=gain, rdNoise=rdnoise,
			                                       retPars=True)
			snr = m['g_psfFlux'][ii] / rmsADU
			fwhm = 2*m['g_fluxRad']*1.1 * 0.455
			skyADUps = skyADU / m['g_expTime'][ii]
			nominal_snr = [ boketc.snr_singleexposure('g',m['psfMag'][i,1],
			                                          m['g_expTime'][i],
			                                          fwhm=fwhm[i],
			                                          skyADU=skyADUps[0],
			                                          profile='gaussian')
			                        for i in ii ]
			nominal_snr = np.array(nominal_snr)
			# revise the ETC calculation using updated gain and RN values,
			# as well as the noise-equivalent-gaussian determined from the
			# pixel area of the PSF
			NEG = np.sqrt(A/(4*np.pi)) * 0.455 * 2.355
			revnominal_snr = [ boketc.snr_singleexposure('g',m['psfMag'][i,1],
			                                             m['g_expTime'][i],
			                                             fwhm=NEG,
			                                             skyADU=skyADUps[0],
			                                             profile='gaussian',
			                                             gain=gain,
			                                             rdnoise=rdnoise)
			                        for i in ii ]
			revnominal_snr = np.array(revnominal_snr)
			objEl = m['g_psfFlux'][ii] * gain 
			est_snr = objEl / np.sqrt(objEl + rmsEl**2)
			sex_snr = m['g_psfFlux'][ii] / m['g_psfFluxErr'][ii]
			ents.extend( [ vals for vals in itertools.izip(ii,
			                                               m['psfMag'][ii,1],
			                                               [A]*len(ii),
			                                               skyADUps,fwhm,
			                                               snr,nominal_snr,
			                                               est_snr,sex_snr,
			                                               revnominal_snr) ] )
		ents = np.array(ents,dtype=[('ii','i4'),('refMag','f4'),
		                            ('psfArea','f4'),('skyADUperSec','f4'),
		                            ('fwhm','f4'),
		                            ('snrRMS','f4'),('snrETC','f4'),
		                            ('snrSky','f4'),('snrSex','f4'),
		                            ('snrETCrev','f4')])
		fitsio.write('cfhtlswide_snr.fits',ents,clobber=(ccdNum==1))

def plot_cfhtls_snr_ratio(snr1='snrRMS',snr2='snrETCrev'):
	hdus = fitsio.FITS('cfhtlswide_snr.fits')
	ccds = [hdu.read() for hdu in hdus[1:]]
	plt.figure()
	for pnum,ccd in enumerate(ccds,start=1):
		ax = plt.subplot(2,2,pnum)
		plt.hexbin(ccd['refMag'],ccd[snr1]/ccd[snr2],
		           extent=(20,23.5,0.5,1.5),cmap=plt.cm.Blues)
		plt.axhline(1,c='r')
		ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
		ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
		ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
		ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
		ax.set_title('CCD%d'%pnum)
	plt.figtext(0.01,0.5,'SNR/SNR(ETC)',va='center',rotation='vertical')
	plt.figtext(0.5,0.01,'g mag (CFHTLS)',ha='center')




##############################################################################
#                                                                            #
#                         Pan-STARRS Medium Deeps                            #
#                                                                            #
##############################################################################

def match_ps1mds(matchRad=2.5):
	raise NotImplementedError
	pstiles = panstarrs_md_tiles(observed=True)
	for field,tiles in pstiles.items():
		stars = fitsio.read(ps1md_starfile(field))
		matches = match_objects(stars,tiles)
		fitsio.write('ps1%s_match.fits'%field,matches,clobber=True)




##############################################################################
#                                                                            #
#                             fake sources                                   #
#                                                                            #
##############################################################################

from astropy.io import fits

def fake_sdss_stars_on_tile(stars,tile,
	                        nresample=200,magrange=(22.0,23.4),
	                        stampSize=25,margin=50,aper='psf',
	                        keepfakes=False,savestars=False):
	pixlo = lambda _x: _x-stampSize/2
	pixhi = lambda _x: _x-stampSize/2 + stampSize
	fakemags = np.zeros(nresample*4,dtype=np.float32)
	fakesnr = -np.ones_like(fakemags)
	if aper=='auto':
		magk,fluxk,errk = 'MAG_AUTO','FLUX_AUTO','FLUXERR_AUTO'
	elif aper=='psf':
		magk,fluxk,errk = 'MAG_PSF','FLUX_PSF','FLUXERR_PSF'
	else:
		raise ValueError
	for ccdNum in range(1,5):
		catpath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
		                       tile['fileName']+'_ccd%d.cat.fits'%ccdNum)
		if not os.path.exists(catpath):
			print ' ... %s does not exist, skipping' % catpath
			continue
		cat = fitsio.read(catpath)
		impath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
		                      tile['fileName']+'_ccd%d.fits'%ccdNum)
		_impath = impath.replace('.fits','_pv.fits')
		fakeim = fits.open(_impath)
		im = fakeim[0].data
		nY,nX = im.shape
		ii = np.where( (stars['ra']>cat['ALPHA_J2000'].min()+3e-3) &
		               (stars['ra']<cat['ALPHA_J2000'].max()-3e-3) &
		               (stars['dec']>cat['DELTA_J2000'].min()+3e-3) &
		               (stars['dec']<cat['DELTA_J2000'].max()-3e-3) )[0]
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
		print 'matched %d/%d stars, max scale factor %.2e' % \
		        (len(m1),len(ii),fscale.max())
		fakex = np.random.randint(margin,nX-margin,nresample)
		fakey = np.random.randint(margin,nY-margin,nresample)
		for x,y,fx,fy,fscl in zip(np.round(cat['X_IMAGE'][m2[jj[rindx]]]),
		                          np.round(cat['Y_IMAGE'][m2[jj[rindx]]]),
		                          fakex,fakey,fscale):
			stamp = im[pixlo(y):pixhi(y),pixlo(x):pixhi(x)]
			im[pixlo(fy):pixhi(fy),pixlo(fx):pixhi(fx)] += fscl*stamp
		fakeimpath = impath.replace('.fits','_fake.fits')
		fakecatpath = fakeimpath.replace('.fits','.cat.fits')
		fakeim.writeto(fakeimpath,clobber=True)
		bokextract.sextract(fakeimpath,frompv=False,redo=True,
		                    withpsf=True,redopsf=False,
		                    psfpath=impath.replace('.fits','.ldac_cat.psf'))
		fakecat = fitsio.read(fakecatpath)
		q1,q2 = srcorXY(fakex,fakey,fakecat['X_IMAGE'],fakecat['Y_IMAGE'],3.0)
		snr = fakecat[fluxk][q2] / fakecat[errk][q2]
		fakemags[nresample*(ccdNum-1):nresample*ccdNum] = fakemag
		fakesnr[nresample*(ccdNum-1):nresample*ccdNum][q1] = snr
		if True:
			zpt = np.median(cat[magk][m2[jj]] - stars['psfMag_g'][ii[m1[jj]]])
			zpt -= 25
			foo = np.where(fakemag[q1] < 22.3)[0]
			offset = np.median((-2.5*np.log10(fakecat[fluxk][q2[foo]]) - zpt) - fakemag[q1[foo]])
			print 'fake star mag offset is ',offset
			fakemags[nresample*(ccdNum-1):nresample*ccdNum] += offset
		if False:
			print ' --------- ZERO POINT CHECK -----------'
			print cat[magk][m2[jj]][:10]
			print -2.5*np.log10(cat[fluxk][m2[jj]])[:10] - zpt
			print stars['psfMag_g'][ii[m1]][:10]
			print ( (-2.5*np.log10(cat[fluxk][m2[jj]])[:10] - zpt) - 
			            stars['psfMag_g'][ii[m1]][:10])
			print -2.5*np.log10(fakecat[fluxk][q2[foo]]) - zpt
			print fakemag[q1[foo]]
			print ( (-2.5*np.log10(fakecat[fluxk][q2[foo]]) - zpt) - 
			         fakemag[q1[foo]] )
			print ( (-2.5*np.log10(fakecat[fluxk][q2[foo]]) - zpt) - 
			         fakemag[q1[foo]] ).mean()
			print snr[foo]
			print 
		if not keepfakes:
			os.unlink(fakeimpath)
			os.unlink(fakecatpath)
		if savestars:
			np.savetxt(fakeimpath.replace('.fits','_stars.dat'),
			   np.vstack([fakemag,fakex,fakey]).transpose(),fmt='%9.3f')
	return fakemags,fakesnr

def fake_ndwfs_stars(grange=(16.0,17.0),**kwargs):
	aper = kwargs.setdefault('aper','psf')
	magrange = kwargs.setdefault('magrange',(22.0,23.4))
	nbins = 5
	medges = np.linspace(magrange[0],magrange[1],nbins+1)
	np.random.seed(1)
	stars = fitsio.read('/global/scratch2/sd/imcgreer/ndwfs/sdss_bootes_gstars.fits')
	fakedir = '/global/scratch2/sd/imcgreer/fakes/'
	stars = stars[(stars['psfMag_g']>grange[0])&(stars['psfMag_g']<grange[1])]
	tiles = ndwfs_tiles(observed=True)
	summaryf = open(fakedir+'fakestars_%s_bytile.dat' % aper,'w')
	summaryf.write('# %4s %1s %8s ' % ('tile','D','utdate'))
	for i in range(nbins):
		summaryf.write('%6.3f ' % ((medges[i]+medges[i+1])/2))
	summaryf.write('\n')
	for ti,tile in enumerate(tiles):
		print 'faking stars in tile %d/%d' % (ti+1,len(tiles))
		mag,snr = fake_sdss_stars_on_tile(stars,tile,**kwargs)
		np.savetxt(fakedir+'fakestars_%s_%05d_%d_%s.dat' % 
		           (aper,tile['tileId'],tile['ditherId'],tile['utDate']),
		           np.vstack([mag,snr]).transpose(),fmt='%8.3f')
		summaryf.write(' %05d %1d %8s ' %
		               (tile['tileId'],tile['ditherId'],tile['utDate']))
		ii = np.digitize(mag,medges)
		# could divide by CCD
		for i in range(nbins):
			jj = np.where(ii==i+1)[0]
			frac = np.sum(snr[jj]>5.0) / float(len(jj))
			summaryf.write('%6.3f ' % frac)
		summaryf.write('\n')
	summaryf.close()



def ndwfs_sdss_matches():
	''' for checking linearity '''
	import basslog
	stars = fitsio.read('/global/scratch2/sd/imcgreer/ndwfs/sdss_bootes_gstars.fits')
	logs = basslog.load_Bok_logs('./logs/')
	tiles = ndwfs_tiles(observed=True)
	tiledb = bass.load_tiledb()
	tid = np.array([int(tid) for tid in tiledb['TID']])
	i1 = 0
	m = np.zeros(1e5,dtype=[('sdss_id','i4'),('sdss_g_mag','f4'),
	                        ('bass_g_mag','f4'),('bass_g_err','f4'),
                            ('bass_expTime','f4'),('bass_skyADU','f4'),
	                        ('bass_airmass','f4'),('bass_ebv','f4'),
	                        ('bass_ccdNum','i4'),('bass_ditherId','i4'),
	                        ('bass_fluxMax','f4'),('bass_FWHM','f4')])
	for ti,tile in enumerate(tiles):
		print 'tile %d/%d [%d]' % (ti+1,len(tiles),i1)
		for ccdNum in range(1,5):
			impath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
			                      tile['fileName']+'_ccd%d.fits'%ccdNum)
			if not os.path.exists(impath):
				print ' ... %s does not exist, skipping' % impath
				continue
			h = fitsio.read_header(impath)
			sky = h['SKYVAL']
			catpath = os.path.join(bass.rdxdir,tile['utDate'],'ccdproc3',
			                       tile['fileName']+'_ccd%d.cat.fits'%ccdNum)
			if not os.path.exists(catpath):
				print ' ... %s does not exist, skipping' % catpath
				continue
			cat = fitsio.read(catpath)
			ii = np.where( (stars['ra']>cat['ALPHA_J2000'].min()+3e-3) &
			               (stars['ra']<cat['ALPHA_J2000'].max()-3e-3) &
			               (stars['dec']>cat['DELTA_J2000'].min()+3e-3) &
			               (stars['dec']<cat['DELTA_J2000'].max()-3e-3) )[0]
			if len(ii)==0:
				print 'no stars found on ccd #',ccdNum
				continue
			m1,m2 = srcor(stars['ra'][ii],stars['dec'][ii],
			              cat['ALPHA_J2000'],cat['DELTA_J2000'],2.5)
			jj = np.where(cat['FLAGS'][m2] == 0)[0]
			i2 = i1 + len(jj)
			m['sdss_id'][i1:i2] = ii[m1[jj]]
			m['sdss_g_mag'][i1:i2] = stars['psfMag_g'][ii[m1[jj]]]
			m['bass_g_mag'][i1:i2] = cat['MAG_PSF'][m2[jj]]
			m['bass_g_err'][i1:i2] = cat['MAGERR_PSF'][m2[jj]]
			m['bass_fluxMax'][i1:i2] = cat['FLUX_MAX'][m2[jj]]
			m['bass_FWHM'][i1:i2] = np.median(cat['FWHM_IMAGE'][m2[jj]])
			m['bass_expTime'][i1:i2] = tile['expTime']
			i = np.where(logs[tile['utDate']]['fileName'] ==
			                tile['fileName'])[0][0]
			m['bass_airmass'][i1:i2] = logs[tile['utDate']]['airmass'][i]
			m['bass_ebv'][i1:i2] = tiledb['EBV'][tid==tile['tileId']][0]
			m['bass_ccdNum'][i1:i2] = ccdNum
			m['bass_ditherId'][i1:i2] = tile['ditherId']
			m['bass_skyADU'][i1:i2] = sky
			i1 = i2
	m = m[:i1]
	outdir = '/project/projectdirs/cosmo/staging/bok/ian/'
	fitsio.write(outdir+'ndwfs_sdss.fits',m,clobber=True)





def get_phototiles_info():
	import boklog
	logs = boklog.load_Bok_logs('./logs/')
	tiledb = bass.load_tiledb()
	tid = np.array([int(tid) for tid in tiledb['TID']])
	ccdNum = 1
	photinfof = open('photo_tiles_info.txt','w')
	photinfof.write('# %6s %10s %7s %7s %7s %10s %8s %7s\n' %
	       ('UTD','file','airmass','E(B-V)','FWHMpix','skyADU','zpt','texp'))
	for ti,tiles in enumerate([cfhtw3_tiles(),ndwfs_tiles()]):
		if ti==0:
			refcat = fitsio.read(cfhtlswide_starfile)
			ii = np.where((refcat['psfMag'][:,1]>17) & 
			              (refcat['psfMag'][:,1]<18.5))[0]
			ref_ra = refcat['ra'][ii]
			ref_dec = refcat['dec'][ii]
			ref_mag = refcat['psfMag'][ii,1]
			#ref_mag = refcat['psfMag'][ii,1] - A_ext['g']*refcat['E(B-V)'][ii]
		else:
			refcat = fitsio.read(bootes_sdss_starfile)
			ii = np.where((refcat['psfMag_g']>16) & 
			              (refcat['psfMag_g']<18.5))[0]
			ref_ra = refcat['ra'][ii]
			ref_dec = refcat['dec'][ii]
			ref_mag = refcat['psfMag_g'][ii]
			#ref_mag = refcat['psfMag_g'][ii] - refcat['extinction_g'][ii]
		for tj,t in enumerate(tiles):
			if t['ditherId'] != 1:
				continue
			# get E(B-V) from tile database
			ebv = tiledb['EBV'][tid==t['tileId']][0]
			# get conditions (airmass,exptime) from observing logs
			try:
				i = np.where(logs[t['utDate']]['fileName']==t['fileName'])[0][0]
			except:
				continue
			airmass = logs[t['utDate']]['airmass'][i]
			exptime = logs[t['utDate']]['expTime'][i]
			# get sky value in ADU from FITS headers
			impath = os.path.join(bass.rdxdir,t['utDate'],'ccdproc3',
			                      t['fileName']+'_ccd%d.fits'%ccdNum)
			h = fitsio.read_header(impath)
			sky = h['SKYVAL']
			# get FWHM and zero point from catalogs
			catpath = os.path.join(bass.rdxdir,t['utDate'],'ccdproc3',
			                       t['fileName']+'_ccd%d.cat.fits'%ccdNum)
			cat = fitsio.read(catpath)
			ii = np.where( (ref_ra>cat['ALPHA_J2000'].min()+3e-3) &
			               (ref_ra<cat['ALPHA_J2000'].max()-3e-3) &
			               (ref_dec>cat['DELTA_J2000'].min()+3e-3) &
			               (ref_dec<cat['DELTA_J2000'].max()-3e-3) )[0]
			if len(ii)==0:
				continue
			m1,m2 = srcor(ref_ra[ii],ref_dec[ii],
			              cat['ALPHA_J2000'],cat['DELTA_J2000'],2)
			if len(m1)==0:
				continue
			m1 = ii[m1]
			ii = np.where(cat['FLAGS'][m2]==0)[0]
			m1,m2 = m1[ii],m2[ii]
			if len(m1)<5:
				continue
			print len(ii),' stars on tile ',t['utDate'],t['fileName']
			fwhm = np.median(cat['FWHM_IMAGE'][m2])
			zpt = 25 - np.median(cat['MAG_AUTO'][m2] - ref_mag[m1]) - \
			         2.5*np.log10(exptime)
			photinfof.write('%8s %10s %7.2f %7.3f %7.2f %10.2f %8.3f %7.1f\n' %
			     (t['utDate'],t['fileName'],airmass,ebv,fwhm,sky,zpt,exptime))
	photinfof.close()

def phototiles_stats(doplots=True):
	import boketc
	gain = boketc.G
	pxscl = boketc.p
	k = boketc.k_ext['g']
	A = boketc.A_ext['g']
	tiledat = ascii_io.read('photo_tiles_info.txt')
	sky_ADUs = tiledat['skyADU'] / tiledat['texp'] 
	sky_eps = sky_ADUs * gain 
	sky_magasec2 = -2.5*np.log10(sky_ADUs*pxscl**-2) + tiledat['zpt']
	print sky_ADUs.mean(),sky_eps.mean(),sky_magasec2.mean()
	zp0 = tiledat['zpt'] - k*(tiledat['airmass']-1) #- A*tiledat['E(B-V)']
	print zp0.mean()
	fwhm_asec = tiledat['FWHMpix'] * pxscl
	if doplots:
		fig = plt.figure(figsize=(8,6))
		ax1 = plt.subplot(2,2,1)
		ax1.hist(zp0)
		#ax1.axvline(boketc.bok_zpt0_am00['g'],c='r',lw=2)
		ax1.axvline(boketc.bok_zpt0_am10['g'],c='r',lw=2)
		ax1 = plt.subplot(2,2,2)
		ax1.hist(sky_magasec2)
		ax1.axvline(boketc.kpno_sky_lun0['g'],c='r',lw=2)
		ax1 = plt.subplot(2,2,3)
		ax1.hist(fwhm_asec)
		ax1.axvline(boketc.bok_medianFWHM['g'],c='r',lw=2)

if __name__=='__main__':
	import sys
	if sys.argv[1]=='match_ndwfs':
		match_ndwfs_stars()
	elif sys.argv[1]=='match_cfhtlswide':
		match_cfhtls_stars(survey='wide')
	elif sys.argv[1]=='fake_ndwfs':
		if len(sys.argv)==2 or 'psf' in sys.argv[2:]:
			aper = 'psf'
		elif 'auto' in sys.argv[2:]:
			aper = 'auto'
		fake_ndwfs_stars(aper=aper)
	elif sys.argv[1]=='photo_info':
		get_phototiles_info()
	else:
		raise ValueError

