#!/usr/bin/env python

import os
import glob
import subprocess

from bass import rdxdir

def sextract(imagepath,frompv=True,withpsf=True,redo=False):
	catpath = imagepath.replace('.fits','.cat.fits')
	if not redo and os.path.exists(catpath):
		print catpath,' exists; skipping'
		return
	_imagepath = imagepath.replace('.fits','_pv.fits') if frompv else imagepath
	if withpsf:
		ldaccatpath = catpath.replace('.cat','.ldac_cat')
		psfpath = ldaccatpath.replace('.fits','.psf')
		cmd = ['sex','-c','config/psfex.sex',
		       '-CATALOG_NAME',ldaccatpath,_imagepath]
		subprocess.call(cmd)
		cmd = ['psfex','-c','config/default.psfex',ldaccatpath]
		subprocess.call(cmd)
		cmd = ['sex','-c','config/default.sex',
		       '-CATALOG_NAME',catpath,
		       '-PSF_NAME',psfpath,_imagepath]
		subprocess.call(cmd)
	else:
		cmd = ['sex','-c','config/default.sex',
		       '-CATALOG_NAME',catpath,_imagepath]
		subprocess.call(cmd)

def sextract_all(**kwargs):
	files = sorted(glob.glob(os.path.join(rdxdir,'2015????','ccdproc3',
	                                      'd????.????_ccd?.wcs')))
	for f in files:
		sextract(f.replace('.wcs','.fits'),**kwargs)


if __name__=='__main__':
	import sys
	redo = 'redo' in sys.argv[1:]
	sextract_all(redo=redo)

