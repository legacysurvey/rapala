#!/usr/bin/env python

import os
import glob
import subprocess

from bass import rdxdir

def sextract(imagepath):
	catpath = imagepath.replace('.fits','.cat.fits')
	if os.path.exists(catpath):
		print catpath,' exists; skipping'
		return
	cmd = ['sex','-c','config/default.sex','-CATALOG_NAME',catpath,imagepath]
	subprocess.call(cmd)

def sextract_all():
	files = sorted(glob.glob(os.path.join(rdxdir,'201504*','ccdproc3',
	                                      'd????.????_ccd?.wcs')))
	for f in files:
		sextract(f.replace('.wcs','.fits'))


if __name__=='__main__':
	sextract_all()

