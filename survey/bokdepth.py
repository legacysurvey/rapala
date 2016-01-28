#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from astropy.stats import sigma_clip
import astropy.io.ascii as ascii_io
import fitsio
try:
	import bass
except ImportError:
	pass

nX,nY = 4032,4096

def make_PSFEx_psf(psfdata,hdr,x_im,y_im):
	assert (hdr['POLNAME1'].strip()=='X_IMAGE' and 
	        hdr['POLNAME2'].strip()=='Y_IMAGE')
	x = (x_im - hdr['POLZERO1'])/hdr['POLSCAL1']
	y = (y_im - hdr['POLZERO2'])/hdr['POLSCAL2']
	deg = hdr['POLDEG1']
	terms = np.array([ y**i * x**j 
	                     for i in range(deg+1) 
	                      for j in range(deg+1) 
	                       if (i+j)<=deg ])
	psfim = psfdata['PSF_MASK'].squeeze() * terms[:,np.newaxis,np.newaxis]
	return psfim.sum(axis=0)

def make_PSFEx_psf_fromfile(psfimf,x_im,y_im):
	'''Given a PSFEx psf file and an (x,y) position on an image, return
	   the corresponding PSF model at that position.
	'''
	psfdata,hdr = fitsio.read(psfimf,header=True)
	return make_PSFEx_psf(psfdata,hdr,x_im,y_im)

def calc_PSF_NEA_grid(psfimf,npts=16):
	'''Given a PSFEx psf file, calculate the Noise Effective Area from the
	   profile weights across a grid of image positions. 'npts' determines
	   the resolution of the grid.
	'''
	nea = np.empty((npts,npts))
	psfdata,hdr = fitsio.read(psfimf,header=True)
	for i in range(npts):
		y = npts/2 + i*(nY//npts)
		for j in range(npts):
			x = npts/2 + j*(nX//npts)
			p = make_PSFEx_psf(psfdata,hdr,x,y)
			p /= p.sum()
			nea[i,j] = np.sum(p**2)**-1
	return nea

def noise_dist_from_psf(imf,psfimf,nsample=1000,_rmhack=False):
	im = fitsio.read(imf)
	psfdata,hdr = fitsio.read(psfimf,header=True)
	if _rmhack:
		# rotated...
		X = 50 + (nY-100) * np.random.rand(nsample)
		Y = 50 + (nX-100) * np.random.rand(nsample)
	else:
		X = 50 + (nX-100) * np.random.rand(nsample)
		Y = 50 + (nY-100) * np.random.rand(nsample)
	psfSum = np.empty(nsample)
	for i in range(nsample):
		psf = make_PSFEx_psf(psfdata,hdr,X[i],Y[i])
		psf /= psf.sum()
		ny,nx = psf.shape
		xi,yi = int(X[i]),int(Y[i])
		i1 = yi - ny//2
		i2 = i1 + ny
		j1 = xi - nx//2
		j2 = j1 + nx
		#stamp = im[max(0,i1):min(nY,i2),max(0,j1):min(nX,j2)]
		stamp = im[i1:i2,j1:j2]
		psfSum[i] = np.sum(psf*stamp)
	return psfSum

def calc_raw_image_background(imagepath,extNum=None,**kwargs):
	margin = kwargs.get('stat_margin',500)
	stat_dtype = [('ampNum','i4'),('skyMean','f4'),('skyMedian','f4'),
	               ('skyVar','f4'),('biasLevel','f4')]
	fits = fitsio.FITS(imagepath)
	if extNum is None:
		extNums = range(1,17)
	elif type(extNum) in [int,str]:
		extNums = [extNum]
	else:
		extNums = extNum
	rv = np.empty(len(extNums),dtype=stat_dtype)
	for i,extn in enumerate(extNums):
		im,bias = improcess(fits[extn],extn,**kwargs)
		pix = sigma_clip(im[margin:-margin,margin:-margin],iters=3)
		extnum = int(fits[extn].get_extname().upper().replace('IM','')) # "IM4"->4
		rv['ampNum'][i] = int(extnum)
		rv['skyMean'][i] = pix.mean()
		rv['skyMedian'][i] = np.ma.median(pix)
		rv['skyVar'][i] = pix.var()
		rv['biasLevel'][i] = bias
	return rv

def calc_sky_all(survey='bass',filt='g'):
	if survey=='bass':
		import basslog
		datadir = os.environ['GSCRATCH']
		master_pixflat = datadir+'/rmreduce/20150315/domeflat_g.fits'
		#master_supersky = datadir+'/rmreduce/20150315/superflt_g.fits'
		logs = basslog.load_Bok_logs('logs/')
		dpfx = ''
	else:
		import boklog
		datadir = os.environ['BOK90PRIMERAWDIR']
		master_pixflat = os.environ['BOK90PRIMEOUTDIR']+'/ut20140315/domeflat_g.fits'
		logs = boklog.load_Bok_logs()
		dpfx = 'ut'
	tile_dtype = [('utDate','S8'),('fileName','S30'),('expTime','f4')]
	tileList = []
	statList = []
	pixflatim = load_flat_im(master_pixflat)
	#superskyim = fitsio.FITS(master_supersky)
	utds = sorted(logs.keys())
	for utd in utds:
		ii = np.where((logs[utd]['imType']=='object') &
		              (logs[utd]['expTime']>30.0) &
		              (logs[utd]['filter']==filt))[0]
		# takes too long to do all of them
		stride = 10 if survey=='bass' else 1
		for _i,i in enumerate(ii[::stride]):
			tile = logs[utd][i]
			print '[%s] tile %d/%d (%s)' % (utd,_i+1,len(ii),tile['fileName'])
			try:
				imagepath = os.path.join(datadir,dpfx+utd,
				                         tile['fileName']+'.fits.gz')
				imstat = calc_raw_image_background(imagepath,
				                                   pixflatim=pixflatim)
				tileList.append((utd,tile['fileName'],tile['expTime']))
				statList.append(imstat)
			except:
				print 'skipping ',tile['fileName']
				continue
	tileData = np.array(tileList,dtype=tile_dtype)
	statData = np.concatenate(statList)
	outfn = 'tile_stats_%s_%s.fits' % (survey,filt)
	fitsio.write(outfn,tileData,clobber=True)
	fitsio.write(outfn,statData)

def calc_depth_tile(tile,**kwargs):
	imagepath = os.path.join(bass.bass_data_dir,
	                         tile['utDate'],tile['fileName']+'.fits.gz')
	imstat = calc_raw_image_background(imagepath,
	                                   extNum=kwargs.get('extNum',4),**kwargs)
	return imstat

def calc_processed_image_rms(imagepath,psfpath,retPars=False,
                             gain=1.375,rdNoise=7.0):
	'''Given a processed image and the PSFEx PSF model, calculate RMS noise
	   level for that image.
	'''
	im,hdr = fitsio.read(imagepath,header=True)
	skyCountsADU = hdr['SKYVAL']
	centerSlice = sigma_clip(im[1024:3072,1008:3024],iters=3,sig=3.5)
	pixRms = centerSlice.std() 
	psfCenter = make_PSFEx_psf_fromfile(psfpath,2048,2016)
	psfCenter /= psfCenter.sum()
	nea = np.sum(psfCenter**2)**-1
	imRmsADU = pixRms * np.sqrt(nea)
	rmsElectrons = np.sqrt( nea * (gain*skyCountsADU + rdNoise**2) )
	if retPars:
		return imRmsADU,rmsElectrons,nea,skyCountsADU
	else:
		return imRmsADU,rmsElectrons

def calc_depth_all():
	obsdb = bass.load_obsdb()
	outf = open('imagestat.dat','w')
	for ti,tile in enumerate(obsdb):
		print 'tile %d/%d (%s)' % (ti+1,len(obsdb),tile['fileName'])
		if tile['filter']=='g':
			try:
				depthstat = calc_depth_tile(tile)
			except:
				print 'skipping ',tile['fileName']
				continue
			outf.write('%04d %8s %10s %6.1f ' % 
			     (ti,tile['utDate'],tile['fileName'],tile['expTime']))
			outf.write('%8.1f %8.1f %8.1f\n' % 
			     (depthstat['bias'],depthstat['medsky'],depthstat['rmssky']))
	outf.close()

def read_image_stats():
	imstat = ascii_io.read('imagestat.dat',
	                       names=['ti','utDate','fileName','expTime',
	                              'biasADU','skyADU','rmsADU'])
	return imstat

def strip_charts():
	imstat = read_image_stats()
	# XXX use per-image zeropoints once available
	zpt = 25.55
	skymag = -2.5*np.log10(imstat['skyADU']/imstat['expTime']*0.455**-2) + zpt
	plt.figure(figsize=(10,4))
	plt.subplots_adjust(0.06,0.13,0.99,0.98)
	ax1 = plt.subplot(111)
	ax1.plot(skymag)
	utds,ii = np.unique(imstat['utDate'],return_index=True)
	for i,utd in zip(ii,utds):
		offset = 1.5 if np.sum(imstat['utDate']==utd) < 30 else 0.2
		ax1.text(i,skymag[i]-offset,imstat['utDate'][i],size=8,
		         rotation='vertical',ha='left',va='bottom')
	ax1.axvspan(ii[utds==20150305],ii[utds==20150313],color='c',alpha=0.5)
	ax1.axhline(22.1,c='r')
	ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
	ax1.set_ylim(22.5,17)
	ax1.set_xlim(-10,6200)
	ax1.set_xlabel('image number')
	ax1.set_ylabel('g sky brightness [mag/asec2]')

def plot_NEA_images(psfimfs,axes=None,npts=16):
	if axes is None:
		plt.figure()
		axes = [plt.subplot(2,2,pnum) for pnum in range(1,5)]
	for ax,ccdNum,psfimf in zip(axes,[1,3,2,4],psfimfs):
		nea = calc_PSF_NEA_grid(psfimf,npts=npts)
		ax.imshow(nea,origin='lower')

if __name__=='__main__':
	#calc_depth_all()
	#calc_sky_all()
	calc_sky_all('sdssrm')

