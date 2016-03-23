#!/usr/bin/env python

import os,shutil
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

def process_images(images,outputDir='./',overwrite=False,cleanup=True):
	tmpFile = 'tmp.fits'
	tmpWcsFile = tmpFile.replace('.fits','_ext1.wcs')
	# read in the reference distortion map
	distortCfg = os.path.join(os.environ['BOKPIPE'],'..','survey',
	                          'config','bokdistort.ahead')
	distortHdrs = read_headers(distortCfg)
	distortWcs1 = WCS(distortHdrs[0])
	x0,y0 = distortWcs1.wcs.crpix  # default tangent plane origin
	xc,yc = 2048,2016              # center of CCD
	# 
	for image in images:
		imFile = os.path.basename(image).replace('.fz','')
		catFile = os.path.join(outputDir,imFile.replace('.fits','.cat.fits'))
		psfFile = os.path.join(outputDir,imFile.replace('.fits','.psf'))
		if os.path.exists(catFile) and not overwrite:
			continue
		combine_ccds([image],output_map=lambda s: tmpFile,
		             clobber=True,_preprocess_function=quick_process_fun)
		if os.path.exists(tmpWcsFile):
			os.remove(tmpWcsFile)
		# use astrometry.net to find the image center in world coords
		solve_bass_image(tmpFile,extns=[1])
		wcsHdr = fits.Header.fromfile(tmpWcsFile)
		wcs1 = WCS(wcsHdr)
		ra0,dec0 = wcs1.all_pix2world(x0,y0,1,ra_dec_order=True)
		ra1,dec1 = wcs1.all_pix2world(xc,yc,1,ra_dec_order=True)
		# set the reference WCS to the updated image origin
		distortWcs1.wcs.crval = [ ra0, dec0 ]
		# find the offset at the center of the image, where the solutions
		#  should be very close
		ra2,dec2 = distortWcs1.all_pix2world(xc,yc,1,ra_dec_order=True)
		# and update the image with these coords, including the offset
		tmpFits = fits.open(tmpFile,mode='update')
		for ccdNum in range(1,5):
			tmpFits[ccdNum].header.update(distortHdrs[ccdNum-1].items())
			tmpFits[ccdNum].header['CRVAL1'] = float(ra0+(ra1-ra2))
			tmpFits[ccdNum].header['CRVAL2'] = float(dec0+(dec1-dec2))
		tmpFits.close()
		# run sextractor+psfex to get object catalogs and PSF models
		sextract(tmpFile,frompv=False,redo=True,
		         withpsf=True,redopsf=True,psfpath=None,onlypsf=True)
		# rename the output files to the original filename
		shutil.move('tmp.ldac_cat.fits',catFile)
		shutil.move('tmp.ldac_cat.psf',psfFile)
		if cleanup:
			os.remove(tmpFile)
			os.remove(tmpWcsFile)
			os.remove('tmp.axy')

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("inputFiles",type=str,nargs='+',
	                    help="input FITS images")
	parser.add_argument("-o","--outputdir",type=str,default='./',
	                    help="output directory")
	parser.add_argument("-R","--redo",action="store_true",
	                    help="reprocess and overwrite existing files")
	parser.add_argument("--noclean",action="store_true",
	                    help="don't delete temporary files")
	args = parser.parse_args()
	process_images(args.inputFiles,
	               outputDir=args.outputdir,
	               overwrite=args.redo,
	               cleanup=not args.noclean)

