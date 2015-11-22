#!/usr/bin/env python

import os
import subprocess

def sextract(imageFile,catFile,psfFile=None,overwrite=False):
	configDir = os.path.join(os.path.split(__file__)[0],'config')
	extractCmd = 'sextractor'
	cmd = [extractCmd,'-c',os.path.join(configDir,'bok_extract.sex')]
	cmd.extend(['-CATALOG_NAME',catFile])
	cmd.extend(['-CATALOG_TYPE','FITS_LDAC'])
	if psfFile is not None:
		raise NotImplementedError
		if psfpath is None:
			ldaccatpath = catpath.replace('.cat','.ldac_cat')
			psfpath = ldaccatpath.replace('.fits','.psf')
		else:
			ldaccatpath = psfpath.replace('.psf','.ldac_cat.fits')
		cmd.extend(['-PSF_NAME',psfpath])
	cmd.extend([
	  '-PARAMETERS_NAME',os.path.join(configDir,'bok_extract.par'),
	  '-FILTER_NAME',os.path.join(configDir,'default.conv'),
	  '-STARNNW_NAME',os.path.join(configDir,'default.nnw'),
	])
	cmd.append(imageFile)
	try:
		subprocess.call(cmd)
	except:
		subprocess.call(['sex']+cmd[1:])
	if psfFile is not None:
		cmd = ['psfex','-c','config/default.psfex',ldaccatpath]
		subprocess.call(cmd)

