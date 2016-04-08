#!/usr/bin/env python

import os,shutil
import glob
from collections import defaultdict
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip
from astropy.table import Table

from bokpipe.bokutil import stats_region
from bokpipe.bokproc import combine_ccds,ampOrder,nominal_gain
from bokpipe.bokoscan import overscan_subtract
from bokpipe.bokastrom import read_headers

from bokastrom import solve_bass_image
from bokextract import sextract
from bass import ampNums,get_amp_index

import ps1cal

gains = { 'IM%d'%ampNum:g for ampNum,g in zip(ampOrder,nominal_gain) }

pixscale = 0.455 # XXX hardcoded

def quick_process_fun(data,hdr,extName):
	im = overscan_subtract(data,hdr,method='mean_value',
	                       reject='sigma_clip',clip_iters=1,
	                       apply_filter=None)
	#im *= gains[extName.upper()]
	return im

def get_ps1_zpastrom(matchedCat,filt,expTime,apNum=1,brightLim=16.0):
	band = filt[-1] # bokr->r
	j = 'gr'.find(band) # index into PS1 arrays
	# read in PS1 color terms
	coeffsf = os.path.join(os.environ['BOKPIPE'],'..','survey',
	                       'config','bok2ps1_%s_coeff.txt'%band)
	ps1colorterms = np.loadtxt(coeffsf)
	# select stars to use for zeropoint measurement
	good = ( (matchedCat['MEDIAN'][:,j] < 25) & 
	         (matchedCat['MEDIAN'][:,j] > brightLim) &
	         (matchedCat['FLAGS'] == 0) )
	# transform PS1 mag to Bok mag using pre-defined color terms
	ps1mag = matchedCat['MEDIAN'][:,j]
	ps1gi = np.diff(matchedCat['MEDIAN'][:,[2,0]],axis=1).squeeze()
	ps1bokmag = ps1mag + np.polyval(ps1colorterms,ps1gi)
	# the intrinsic Bok magnitude
	bokflux = np.ma.array(matchedCat['FLUX_APER'][:,apNum],
	                      mask=((matchedCat['FLUX_APER'][:,apNum]<=0)|~good))
	bokmag = -2.5*np.ma.log10(bokflux/expTime)
	# find the average zeropoint over the full image
	dmag = sigma_clip(ps1bokmag-bokmag)
	zpIm = dmag.mean()
	nim = (~dmag.mask).sum()
	#
	cosdec = np.cos(np.radians(matchedCat['DEC']))
	dra = (matchedCat['RA'] - matchedCat['ALPHA_J2000'])*cosdec
	ddec = matchedCat['DEC'] - matchedCat['DELTA_J2000']
	# now find average zeropoints over each CCD and each amp
	zpCCD = np.zeros(4,dtype=np.float32)
	zpAmp = np.zeros((4,4),dtype=np.float32)
	raoff = np.zeros(4,dtype=np.float32)
	decoff = np.zeros(4,dtype=np.float32)
	seeing = np.zeros(4,dtype=np.float32)
	medgicolor = np.zeros(4,dtype=np.float32)
	nccd = np.zeros(4,dtype=np.int32)
	namp = np.zeros((4,4),dtype=np.int32)
	ampIndex = get_amp_index(matchedCat['X_IMAGE'],matchedCat['Y_IMAGE'])
	for ccdNum in range(1,5):
		ii = np.where(matchedCat['ccdNum']==ccdNum)[0]
		dmag = sigma_clip(ps1bokmag[ii]-bokmag[ii])
		zpCCD[ccdNum-1] = dmag.mean()
		nccd[ccdNum-1] = (~dmag.mask).sum()
		raoff[ccdNum-1] = np.median(dra[ii[~dmag.mask]])*3600
		decoff[ccdNum-1] = np.median(ddec[ii[~dmag.mask]])*3600
		seeing = np.median(matchedCat['FWHM_IMAGE'][ii[~dmag.mask]])*pixscale
		medgicolor[ccdNum-1] = np.median(ps1gi[ii[~dmag.mask]])
		for ai in range(4):
			ii = np.where((matchedCat['ccdNum']==ccdNum) & 
			              (ampIndex==ai))[0]
			dmag = sigma_clip(ps1bokmag[ii]-bokmag[ii])
			zpAmp[ccdNum-1,ai] = dmag.mean()
			namp[ccdNum-1,ai] = (~dmag.mask).sum()
	return dict(zpim=zpIm,zpccd=zpCCD,zpamp=zpAmp,
	            raoff=raoff,decoff=decoff,
	            seeing=seeing,medgicolor=medgicolor,
	            nstarsim=nim,nstarsccd=nccd,nstarsamp=namp)

def _hextract_wcs(hdus):
	tmp = defaultdict(list)
	wcskeys = ['CRPIX1','CRPIX2','CRVAL1','CRVAL2',
	           'CD1_1','CD1_2','CD2_1','CD2_2'] 
	for ccdNum in range(1,5):
		for k in wcskeys:
			tmp[k.lower()].append(hdus[ccdNum].header[k])
	return {k:np.array(v) for k,v in tmp.items()}

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
	ccdcenterpix = stats_region('ccd_central_quadrant')
	# 
	for image in images:
		imFile = os.path.basename(image).replace('.fz','')
		catFile = os.path.join(outputDir,imFile.replace('.fits','.cat.fits'))
		psfFile = os.path.join(outputDir,imFile.replace('.fits','.psf'))
		metaFile = os.path.join(outputDir,imFile.replace('.fits','.meta'))
		# check if all the output files exists 
		if ( os.path.exists(catFile) and os.path.exists(psfFile) and
		       os.path.exists(metaFile) and not overwrite ):
			continue
		combine_ccds([image],output_map=lambda s: tmpFile,
		             clobber=True,_preprocess_function=quick_process_fun)
		if os.path.exists(tmpWcsFile):
			os.remove(tmpWcsFile)
		# store the filter and exptime now before tmp file gets deleted
		hdr0 = fits.getheader(tmpFile,0)
		filt = hdr0['FILTER']
		expTime = hdr0['EXPTIME']
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
		avsky = np.zeros(4,dtype=np.float32)
		#avskyamp = np.zeros((4,4),dtype=np.float32)
		for ccdNum in range(1,5):
			tmpFits[ccdNum].header.update(distortHdrs[ccdNum-1].items())
			tmpFits[ccdNum].header['CRVAL1'] = float(ra0+(ra1-ra2))
			tmpFits[ccdNum].header['CRVAL2'] = float(dec0+(dec1-dec2))
			skypix = tmpFits[ccdNum].data[ccdcenterpix]
			avsky[ccdNum-1] = sigma_clip(skypix).mean()
		hdrWcs = _hextract_wcs(tmpFits)
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
		# calculate the zeropoint from PS1
		ps1m = ps1cal.match_ps1(catFile,isldac=True)
		zps = get_ps1_zpastrom(ps1m,filt,expTime)
		# join all the meta-data
		metadata = dict(zps.items()+hdrWcs.items())
		metadata['avsky'] = avsky
		np.savez(metaFile,**metadata)

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

