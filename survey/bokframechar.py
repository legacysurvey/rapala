#!/usr/bin/env python

import os
import glob

from astropy.io import fits
from astropy.wcs import WCS

from bokpipe.bokproc import combine_ccds,ampOrder,nominal_gain
from bokpipe.bokoscan import overscan_subtract
from bokpipe.bokastrom import read_headers

from bokastrom import solve_bass_image
from bokextract import sextract

gains = { 'IM%d'%ampNum:g for ampNum,g in zip(ampOrder,nominal_gain) }

def quick_process_fun(data,hdr,extName):
	im = overscan_subtract(data,hdr,method='mean_value',
	                       reject='sigma_clip',clip_iters=1,
	                       apply_filter=None)
	im *= gains[extName.upper()]
	return im

def process_images(images):
	tmpFile = 'tmp.fits'
	tmpHdrFile = 'tmp.ahead'
	tmpWcsFile = tmpFile.replace('.fits','_ext1.wcs')
	distortCfg = os.path.join(os.environ['BOKPIPE'],'..','survey',
	                          'config','bokdistort.ahead')
	distortHdrs = read_headers(distortCfg)
	for image in images:
		combine_ccds([image],output_map=lambda s: tmpFile,
		             clobber=True,_preprocess_function=quick_process_fun)
		if os.path.exists(tmpWcsFile):
			os.remove(tmpWcsFile)
		solve_bass_image(tmpFile,extns=[1])
		wcsHdr = fits.Header.fromfile(tmpWcsFile)
		wcsFit = WCS(wcsHdr)
		for i,h in enumerate(distortHdrs):
			if i==0:
				w = WCS(h)
				ra1,dec1 = w.all_pix2world(0,0,0,ra_dec_order=True)
				ra2,dec2 = wcsFit.all_pix2world(0,0,0,ra_dec_order=True)
				dra = ra2 - ra1
				ddec = dec2 - dec1
			h['CRVAL1'] += dra
			h['CRVAL2'] += ddec
			h.totextfile('tmp%d.ahead' % (i+1),clobber=True,endcard=True)
		with open(tmpHdrFile,'w') as outf:
			for hdrfn in ['tmp%d.ahead'%n for n in range(1,5)]:
				with open(hdrfn) as hdrf:
					outf.write(hdrf.read())
		sextract(tmpFile,frompv=False,redo=True,
		         withpsf=True,redopsf=True,psfpath=None,onlypsf=True)

if __name__=='__main__':
	import sys
	process_images(sys.argv[1:])

