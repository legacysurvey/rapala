#!/usr/bin/env python

import os
import shutil
import subprocess
import numpy as np
from astropy.table import Table,vstack
from astropy.wcs import WCS
import fitsio

from .bokutil import load_mask
from .bokastrom import read_headers

try:
	import sep
except ImportError:
	pass

configDir = os.path.join(os.path.split(__file__)[0],'config')

def sextract(imageFile,catFile,psfFile=None,clobber=False,full=False,
             verbose=0,**kwargs):
	if not clobber and os.path.exists(catFile):
		if verbose > 0:
			print catFile,' already exists, skipping'
		return
	cmd = ['sex','-c',os.path.join(configDir,'bok_extract.sex')]
	cmd.extend(['-CATALOG_NAME',catFile])
	if full:
		cmd.extend(['-PARAMETERS_NAME',
		            os.path.join(configDir,'bok_catalog.par')])
		cmd.extend(['-PHOT_APERTURES',"7.5,15.0,22.0"])
		cmd.extend(['-CATALOG_TYPE','FITS_1.0'])
		if psfFile is None:
			psfFile = catFile.replace('.fits','.psf')
		cmd.extend(['-PSF_NAME',psfFile])
	else:
		cmd.extend(['-PARAMETERS_NAME',
		            os.path.join(configDir,'bok_extract.par')])
		cmd.extend(['-CATALOG_TYPE','FITS_LDAC'])
	cmd.extend([
	  '-FILTER_NAME',os.path.join(configDir,'default.conv'),
	  '-STARNNW_NAME',os.path.join(configDir,'default.nnw'),
	])
	if verbose >= 5:
		cmd.extend(['-VERBOSE_TYPE','FULL'])
	elif verbose >= 2:
		cmd.extend(['-VERBOSE_TYPE','NORMAL'])
	elif verbose > 0:
		print 'generating wcs catalog for ',imageFile
	if verbose >= 2:
		outdev = None
	else:
		outdev = open(os.devnull,'w')
	cmd.append(imageFile)
	if verbose > 1:
		print ' '.join(cmd)
	subprocess.call(cmd,stdout=outdev)

def run_psfex(catFile,psfFile=None,clobber=False,verbose=0,**kwargs):
	defPsfFile = catFile.replace('.fits','.psf')
	if psfFile is not None and psfFile != defPsfFile:
		rename = True
	else:
		if psfFile is None:
			psfFile = defPsfFile
		rename = False
	if not clobber and os.path.exists(psfFile):
		if verbose > 0:
			print psfFile,' already exists, skipping'
		return
	cmd = ['psfex','-c',os.path.join(configDir,'bok90.psfex'),catFile]
	if verbose > 2:
		cmd.extend(['-VERBOSE_TYPE','FULL'])
	elif verbose > 1:
		cmd.extend(['-VERBOSE_TYPE','NORMAL'])
	elif verbose > 0:
		print 'generating psf for ',catFile
	if verbose > 1:
		print ' '.join(cmd)
	subprocess.call(cmd)
	if rename:
		shutil.move(defPsfFile,psfFile)

def aper_phot(image,hdr,ra,dec,aperRad,mask=None,varim=None,edge_buf=5,
              **kwargs):
	bkgmethod = kwargs.get('background')
	satVal = hdr['SATUR']
	msk = mask | (image > satVal)
	w = WCS(hdr)
	foot = w.calc_footprint()
	raMin,raMax = foot[:,0].min(),foot[:,0].max()
	decMin,decMax = foot[:,1].min(),foot[:,1].max()
	ii = np.where((ra>raMin)&(ra<raMax)&(dec>decMin)&(dec<decMax))[0]
	x,y = w.wcs_world2pix(ra[ii],dec[ii],0)
	ii2 = np.where((x>edge_buf) & (y>edge_buf) & 
	               (x<image.shape[1]-edge_buf) & 
	               (y<image.shape[0]-edge_buf))[0]
	x = x[ii2]
	y = y[ii2]
	ii = ii[ii2]
	nObj,nAper = len(ii),len(aperRad)
	cts = np.empty((nObj,nAper),dtype=np.float32)
	ctserr = np.empty((nObj,nAper),dtype=np.float32)
	flags = np.empty((nObj,nAper),dtype=np.int32)
	nmasked = np.empty(nObj,dtype=np.int32)
	pkvals = np.empty(nObj,dtype=np.float32)
	# compute global background fit if necessary
	if varim is None or bkgmethod=='global':
		bkg = sep.Background(image)
	if varim is None:
		varim = bkg.globalrms**2
	im = image
	skyann = None
	if bkgmethod == 'global':
		im = image - bkg.back()
	elif bkgmethod == 'local':
		skyann = (25.,40.)
	# perform aperture photometry
	for j,aper in enumerate(aperRad):
		c,crms,f = sep.sum_circle(im,x,y,aper,mask=msk,
		                          gain=1.0,bkgann=skyann,var=varim)
		cts[:,j] = c
		ctserr[:,j] = crms
		flags[:,j] = f
	# clean up bad values
	bad = np.isnan(cts) | np.isnan(ctserr)
	cts[bad] = 0
	ctserr[bad] = 0
	flags[bad] |= 16
	# additional masking and diagnostics
	for i,(_x,_y) in enumerate(zip(x.astype(int),y.astype(int))):
		nmasked[i] = msk[_y-3:_y+3+1,_x-3:_x+3+1].sum() 
		pkvals[i] = image[_y-3:_y+3+1,_x-3:_x+3+1].max() 
		if pkvals[i] > satVal:
			flags[i] |= 4
	rv = (x,y,ii,cts,ctserr,flags,nmasked,pkvals)
	return rv

def aper_phot_image(imageFile,ra,dec,aperRad,badPixMask=None,varIm=None,
	                aHeadFile=None,**kwargs):
	from astropy.io.fits import getheader
	maskIsWhtMap = kwargs.get('mask_is_weight_map',True)
	maskType = kwargs.pop('mask_type','gtzero')
	if maskType is None and maskIsWhtMap:
		maskType = 'zerobad'
	phot_cols = ('x','y','objId','counts','countsErr',
	             'flags','nMasked','peakCounts')
	phot_dtype = ('f4','f4','i4','f4','f4','i4','i4','f4')
	phot_cols += ('ccdNum',)
	phot_dtype += ('i4',)
	if aHeadFile is not None:
		wcsHdrs = read_headers(aHeadFile)
	tabs = []
	fitsData = fitsio.FITS(imageFile)
	for i,hdu in enumerate(fitsData[1:]):
		im = hdu.read()
		extn = hdu.get_extname()
		# argh... header need to be from astropy, fix this!!!
		hdr = getheader(imageFile,extn)
		if aHeadFile is not None:
			for k,v in wcsHdrs[i].items():
				hdr[k] = v
		if badPixMask is None:
			mask = None
		else:
			mask = load_mask(badPixMask[extn],maskType)
		varim = varIm[extn] if varIm is not None else None
		phot = aper_phot(im,hdr,ra,dec,aperRad,mask=mask,varim=varim,**kwargs)
		n = len(phot[0])
		if n==0:
			continue
		# add ccd number
		phot += (np.repeat(i+1,n),)
		t = Table(phot,names=phot_cols,dtype=phot_dtype)
		tabs.append(t)
	if len(tabs)==0:
		return None
	else:
		phot = vstack(tabs)
		return phot

