#!/usr/bin/env python

import os
import glob
import subprocess

def sextract(imagepath):
	catpath = imagepath.replace('.fits','.cat.fits')
	if os.path.exists(catpath):
		print catpath,' exists; skipping'
		return
	cmd = ['sex','-c','config/default.sex','-CATALOG_NAME',catpath,imagepath]
	subprocess.call(cmd)

def sextract_all():
	rdxdir = os.path.join(os.environ['GSCRATCH'],'rmreduce')
	files = sorted(glob.glob(os.path.join(rdxdir,'2015*','ccdproc3',
	                                      'd????.????_ccd?.wcs')))
	for f in files:
		sextract(f.replace('.wcs','.fits'))


if __name__=='__main__':
	sextract_all()

