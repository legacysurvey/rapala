#!/usr/bin/env python

import os
import shutil
import subprocess
import numpy as np
from astropy.table import Table,vstack
import fitsio

from .bokastrom import read_headers,wcs_from_header

try:
	import sep
except ImportError:
	print 'sep is not installed, cannot perform aperture photometry routines'

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
	if verbose > 2:
		cmd.extend(['-VERBOSE_TYPE','FULL'])
	elif verbose > 1:
		cmd.extend(['-VERBOSE_TYPE','NORMAL'])
	elif verbose > 0:
		print 'generating wcs catalog for ',imageFile
	cmd.append(imageFile)
	if verbose > 1:
		print ' '.join(cmd)
	subprocess.call(cmd)

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

def aper_phot(image,hdr,ra,dec,aperRad,badPixMask,edge_buf=5,**kwargs):
	w = wcs_from_header(hdr)
	foot = w.calc_footprint()
	raMin,raMax = foot[:,0].min(),foot[:,0].max()
	decMin,decMax = foot[:,1].min(),foot[:,1].max()
	ii = np.where((ra>raMin)&(ra<raMax)&(dec>decMin)&(dec<decMax))[0]
	x,y = w.wcs_world2pix(ra[ii],dec[ii],0)
	ii2 = np.where((x>edge_buf) & (y>edge_buf) & 
	               (x<4096-edge_buf) & (y<4032-edge_buf))[0]
	x = x[ii2]
	y = y[ii2]
	ii = ii[ii2]
	#bkg = sep.Background(image)
	nObj,nAper = len(ii),len(aperRad)
	cts = np.empty((nObj,nAper),dtype=np.float32)
	ctserr = np.empty((nObj,nAper),dtype=np.float32)
	flags = np.empty((nObj,nAper),dtype=np.int32)
	# only data type sep appears to accept
	_mask = badPixMask.astype(np.int32)
	for j,aper in enumerate(aperRad):
		c,cvar,f = sep.sum_circle(image,x,y,aper,mask=_mask,
		                          gain=1.0,bkgann=(25.,35.))
		cts[:,j] = c
		ctserr[:,j] = np.sqrt(cvar)
		flags[:,j] = f
	return x,y,ii,cts,ctserr,flags

def aper_phot_image(imageFile,ra,dec,aperRad,badPixMask,
	                aHeadFile=None,**kwargs):
	from astropy.io.fits import getheader
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
		mask = badPixMask[extn].read().astype(np.bool)
		phot = aper_phot(im,hdr,ra,dec,aperRad,mask,**kwargs)
		n = len(phot[0])
		if n==0:
			continue
		# add ccd number
		phot += (np.repeat(i+1,n),)
		t = Table(phot,
		        names=('x','y','objId','counts','countsErr','flags','ccdNum'),
		        dtype=('f4','f4','i4','f4','f4','i4','i4'))
		tabs.append(t)
	if len(tabs)==0:
		return None
	else:
		phot = vstack(tabs)
		return phot
