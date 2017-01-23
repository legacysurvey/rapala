#!/usr/bin/env python

import os,shutil
import glob
from collections import defaultdict
import numpy as np
import tempfile
import multiprocessing
from functools import partial

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
	                       'config','bok2ps1_%s_gicoeff.dat'%band)
	ps1colorterms = np.loadtxt(coeffsf)
	# select stars to use for zeropoint measurement
	good = ( (matchedCat['MEDIAN'][:,j] < 19) & 
	         (matchedCat['MEDIAN'][:,j] > brightLim) &
	         (matchedCat['FLAGS'] == 0) )
	# transform PS1 mag to Bok mag using pre-defined color terms
	ps1mag = matchedCat['MEDIAN'][:,j]
	ps1gi = np.diff(matchedCat['MEDIAN'][:,[2,0]],axis=1).squeeze()
	ps1bokmag = ps1mag + np.polyval(ps1colorterms,ps1gi)
	# the intrinsic Bok magnitude
	flux = matchedCat['FLUX_APER']
	if len(flux.shape) > 1:
		flux = flux[:,apNum]
	bokflux = np.ma.array(flux,mask=((flux<=0)|~good))
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
	zprmsCCD = np.zeros(4,dtype=np.float32)
	zpAmp = np.zeros((4,4),dtype=np.float32)
	raoff = np.zeros(4,dtype=np.float32)
	decoff = np.zeros(4,dtype=np.float32)
	seeing = np.zeros(4,dtype=np.float32)
	avsky = np.zeros(4,dtype=np.float32)
	medgicolor = np.zeros(4,dtype=np.float32)
	nccd = np.zeros(4,dtype=np.int32)
	namp = np.zeros((4,4),dtype=np.int32)
	ampIndex = get_amp_index(matchedCat['X_IMAGE'],matchedCat['Y_IMAGE'])
	for ccdNum in range(1,5):
		ii = np.where(matchedCat['ccdNum']==ccdNum)[0]
		if len(ii)<10:
			continue
		dmag = sigma_clip(ps1bokmag[ii]-bokmag[ii])
		zpCCD[ccdNum-1] = dmag.mean()
		zprmsCCD[ccdNum-1] = dmag.std()
		nccd[ccdNum-1] = (~dmag.mask).sum()
		raoff[ccdNum-1] = np.median(dra[ii[~dmag.mask]])*3600
		decoff[ccdNum-1] = np.median(ddec[ii[~dmag.mask]])*3600
		seeing[ccdNum-1] = np.median(matchedCat['FWHM_IMAGE'][ii[~dmag.mask]])*pixscale
		medgicolor[ccdNum-1] = np.median(ps1gi[ii[~dmag.mask]])
		avsky[ccdNum-1] = np.median(matchedCat['BACKGROUND'][ii])
		for ai in range(4):
			ii = np.where((matchedCat['ccdNum']==ccdNum) & 
			              (ampIndex==ai))[0]
			if len(ii) > 5:
				dmag = sigma_clip(ps1bokmag[ii]-bokmag[ii])
				zpAmp[ccdNum-1,ai] = dmag.mean()
				namp[ccdNum-1,ai] = (~dmag.mask).sum()
	return dict(zpim=zpIm,zpccd=zpCCD,zpamp=zpAmp,zprmsCCD=zprmsCCD,
	            raoff=raoff,decoff=decoff,
	            seeing=seeing,medgicolor=medgicolor,avsky=avsky,
	            nstarsim=nim,nstarsccd=nccd,nstarsamp=namp)

def _hextract_wcs(hdus,fromldac=False):
	tmp = defaultdict(list)
	wcskeys = ['CRPIX1','CRPIX2','CRVAL1','CRVAL2',
	           'CD1_1','CD1_2','CD2_1','CD2_2'] 
	for ccdNum in range(1,5):
		if fromldac:
			hdrlist = hdus[1+2*(ccdNum-1)].data[0][0]
			hdrstr = ''.join([c.ljust(80,' ') for c in hdrlist])
			hdr = fits.Header.fromstring(hdrstr)
		else:
			hdr = hdus[ccdNum].header
		for k in wcskeys:
			tmp[k.lower()].append(hdr[k])
	return {k:np.array(v) for k,v in tmp.items()}

def process_raw_images(images,outputDir='./',overwrite=False,cleanup=True,
	                   searchrad=0.25,ignore_missing_wcs=True):
	_tmpf = tempfile.NamedTemporaryFile()
	tmpFile = _tmpf.name+'.fits'
	#
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
		print 'processing ',image
		imFile = os.path.basename(image).replace('.fz','')
		wcsFile = os.path.join(outputDir,imFile.replace('.fits','_ext1.wcs'))
		catFile = os.path.join(outputDir,imFile.replace('.fits','.cat.fits'))
		psfFile = os.path.join(outputDir,imFile.replace('.fits','.psf'))
		metaFile = os.path.join(outputDir,imFile.replace('.fits','.meta'))
		if os.path.exists(tmpWcsFile):
			os.remove(tmpWcsFile)
		# check if the combined image is needed
		newimg = False
		filt = None
		if ( (not os.path.exists(wcsFile) and not ignore_missing_wcs)
		      or not os.path.exists(catFile) 
		      or not os.path.exists(psfFile) 
		     or overwrite ):
			try:
				combine_ccds([image],output_map=lambda s: tmpFile,
				             clobber=True,
				             _preprocess_function=quick_process_fun)
				newimg = True
			except IOError:
				print image,' COMBINE FAILED!!!'
				continue
			else:
				# crashes on cori/edison in memmap unless it is off
				tmpFits = fits.open(tmpFile,mode='update',memmap=False)
				filt = tmpFits[0].header['FILTER']
				expTime = tmpFits[0].header['EXPTIME']
		if filt is None:
			tmphdr = fits.getheader(image,0)
			filt = tmphdr['FILTER']
			expTime = tmphdr['EXPTIME']
		wcsHdr = None
		hdrWcs = None
		if not os.path.exists(wcsFile) and newimg:
			# offset from the focal plane center to the CCD center
			dec = tmpFits[1].header['CRVAL2'] + (59+yc)*0.455/3600
			cosdec = np.cos(np.radians(dec))
			ra = tmpFits[1].header['CRVAL1'] + (182+xc)*0.455*cosdec/3600
			# use astrometry.net to find the image center in world coords
			solve_bass_image(tmpFile,extns=[1],ra=ra,dec=dec,
			                 searchrad=searchrad)
#			# trying here to send in sextractor catalog instead of finding srcs
#			sextract(tmpFile,frompv=False,redo=True,
#			         withpsf=True,redopsf=True,psfpath=None,onlypsf=True)
#			solve_bass_image('tmp.ldac_cat.fits',extns=2,ra=ra,dec=dec,
#			                 fromcat=True)
			try:
				wcsHdr = fits.Header.fromfile(tmpWcsFile)
			except IOError:
				print image,' WCS FAILED!!!'
				continue
		if newimg:
			if not wcsHdr:
				# pull in the previous WCS solution
				wcsHdr = fits.Header.fromfile(wcsFile)
			wcs1 = WCS(wcsHdr)
			ra0,dec0 = wcs1.all_pix2world(x0,y0,1,ra_dec_order=True)
			ra1,dec1 = wcs1.all_pix2world(xc,yc,1,ra_dec_order=True)
			# set the reference WCS to the updated image origin
			distortWcs1.wcs.crval = [ ra0, dec0 ]
			# find the offset at the center of the image, where the solutions
			#  should be very close
			ra2,dec2 = distortWcs1.all_pix2world(xc,yc,1,ra_dec_order=True)
			# and update the image with these coords, including the offset
			avrawgain = np.zeros(4,dtype=np.float32)
			avsky = np.zeros(4,dtype=np.float32)
			#avskyamp = np.zeros((4,4),dtype=np.float32)
			for ccdNum in range(1,5):
				tmpFits[ccdNum].header.update(distortHdrs[ccdNum-1].items())
				tmpFits[ccdNum].header['CRVAL1'] = float(ra0+(ra1-ra2))
				tmpFits[ccdNum].header['CRVAL2'] = float(dec0+(dec1-dec2))
				skypix = tmpFits[ccdNum].data[ccdcenterpix]
				avsky[ccdNum-1] = float(sigma_clip(skypix,iters=1).mean())
			hdrWcs = _hextract_wcs(tmpFits)
			tmpFits.close()
			if os.path.exists(tmpWcsFile):
				shutil.move(tmpWcsFile,wcsFile)
			# ... at this point, if we need a new image, it's because either
			#     the catalog or psf files are missing, or both. the sextract()
			#     call generates both so just run it here
			# run sextractor+psfex to get object catalogs and PSF models
			sextract(tmpFile,frompv=False,redo=True,
			         withpsf=True,redopsf=True,psfpath=None,onlypsf=True)
			# rename the output files to the original filename
			shutil.move(tmpFile.replace('.fits','.ldac_cat.fits'),catFile)
			shutil.move(tmpFile.replace('.fits','.ldac_cat.psf'),psfFile)
		else:
			# at this point it means we're only missing the zeropoints and
			# want to recalculate them. but we need to fill in the WCS and
			# other info that would have been found above. can obtain it from
			# the information stored in the sextractor catalogs.
			catFits = fits.open(catFile)
			hdrWcs = _hextract_wcs(catFits,fromldac=True)
			avrawgain = np.zeros(4,dtype=np.float32)
			# this is kind of hacky (but then so is everything else), but
			# since we haven't processed the image, rather than calculate a
			# rough sky estimate from pixels, use the background estimation
			# from sextractor
			avsky = np.array([np.median(catFits[n].data['BACKGROUND'])
			                    for n in range(2,9,2)]).astype(np.float32)
		if not os.path.exists(metaFile):
			# calculate the zeropoint from PS1
			ps1m = ps1cal.match_ps1(catFile,isldac=True)
			if ps1m is None:
				print image,' PS1CAL FAILED!!!'
				continue
			zps = get_ps1_zpastrom(ps1m,filt,expTime)
			# join all the meta-data
			metadata = dict(zps.items()+hdrWcs.items())
			metadata['arawgain'] = avrawgain
			metadata['avsky'] = avsky
			np.savez(metaFile,**metadata)
	if cleanup:
		for f in [tmpFile,tmpWcsFile,tmpFile.replace('.fits','.axy')]:
			try:
				os.remove(f)
			except:
				pass
	_tmpf.close()

def _process_image(image,outputDir='./',overwrite=False,naocver=False,
                   redoifnocal=True):
	print 'processing ',image
	imDir,imFile = os.path.split(image)
	imDir,utDir = os.path.split(imDir)
	outputDir = os.path.join(outputDir,utDir)
	isfz = imFile.endswith('.fz')
	imFile = imFile.replace('.fz','')
	catFile = os.path.join(outputDir,imFile.replace('.fits','.cat.fits'))
	psfFile = os.path.join(outputDir,imFile.replace('.fits','.psf'))
	ps1File = os.path.join(outputDir,imFile.replace('.fits','.ps1match.fits'))
	if ( not os.path.exists(catFile) or not os.path.exists(psfFile) 
	      or overwrite
	      or (not os.path.exists(ps1File) and redoifnocal) ):
		if isfz:
			_image = os.path.join(outputDir,imFile)
			print 'funpack %s -O %s' % (image,_image)
			os.system('funpack -O %s %s' % (_image,image))
			image = _image
		sextract(image,frompv=False,redo=True,
		         withpsf=True,redopsf=True,psfpath=psfFile,onlypsf=True)
		shutil.move(psfFile.replace('.psf','.ldac_cat.fits'),catFile)
		shutil.move(psfFile.replace('.psf','.ldac_cat.psf'),psfFile)
		if isfz:
			os.remove(image)
	if not os.path.exists(ps1File) or overwrite:
		ps1m = ps1cal.match_ps1(catFile,isldac=True,singleccd=naocver)
		if ps1m is None:
			print image,' PS1CAL FAILED!!!'
		else:
			if naocver:
				ps1m['ccdNum'] = int(imFile[-6])
			ps1m.write(ps1File,overwrite=True)

def process_image(image,**kwargs):
	try:
		_process_image(image,**kwargs)
	except:
		print image,' FAILED TO PROCESS'

def process_images(images,nproc=1,**kwargs):
	# have to make the dirs before splitting into subprocesses
	outputDir = kwargs.get('outputDir','.')
	dirs = { os.path.join(outputDir,os.path.split(os.path.split(f)[0])[1]):0
	           for f in images }
	for d in dirs.keys():
		if not os.path.exists(d):
			os.makedirs(d)
	if True:
		# prune all frames that already completed to the PS1 matching step
		allps1m = glob.glob(os.path.join(outputDir,'*','*ps1match.*'))
		allps1m = [ os.path.basename(f).replace('.ps1match','') 
		                   for f in allps1m ]
		images = [ im for im in images 
		               if os.path.basename(im) not in allps1m ]
	print len(images), ' frames to process'
	pool = multiprocessing.Pool(nproc)
	_proc = partial(process_image,**kwargs)
	pool.map(_proc,images)

def process_idm_images(images,outputDir='./',overwrite=False):
	# 
	for image in images:
		print 'processing ',image
		imFile = os.path.basename(image).replace('.fz','')
		catFile = image.replace('.fits','.cat.fits')
		psfFile = image.replace('.fits','.psf.fits')
		metaFile = os.path.join(outputDir,imFile.replace('.fits','.meta'))
		# check if all the output files exists 
		if ( os.path.exists(catFile) and os.path.exists(psfFile) and
		       os.path.exists(metaFile) and not overwrite ):
			continue
		tmpFits = fits.open(image)
		hdrWcs = _hextract_wcs(tmpFits)
		filt = tmpFits[0].header['FILTER']
		expTime = tmpFits[0].header['EXPTIME']
		avrawgain = np.zeros(4,dtype=np.float32)
		avsky = np.zeros(4,dtype=np.float32)
		#avskyamp = np.zeros((4,4),dtype=np.float32)
		for ccdNum in range(1,5):
			hdr = tmpFits[ccdNum].header
			# include all of the gain corrections applied
			gains = [(hdr['GAIN%02dA'%ampNum]) 
			           * (hdr['GAIN%02dB'%ampNum]) 
			            * hdr['CCDGAIN'] 
			             for ampNum in range(4*(ccdNum-1)+1,4*ccdNum+1)]
			avrawgain[ccdNum-1] = np.mean(gains)
			skyvals = [(hdr['SKY%02dB'%ampNum]) 
			             for ampNum in range(4*(ccdNum-1)+1,4*ccdNum+1)]
			avsky[ccdNum-1] = np.mean(skyvals) / np.mean(gains)
		# calculate the zeropoint from PS1
		ps1m = ps1cal.match_ps1(catFile,isldac=False)
		zps = get_ps1_zpastrom(ps1m,filt,expTime)
		# join all the meta-data
		metadata = dict(zps.items()+hdrWcs.items())
		metadata['arawgain'] = avrawgain
		metadata['avsky'] = avsky
		np.savez(metaFile,**metadata)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("inputFiles",type=str,nargs='+',
	                    help="input FITS images")
	parser.add_argument("-o","--outputdir",type=str,default='./',
	                    help="output directory")
	parser.add_argument("-r","--raw",action="store_true",
	                    help="input is raw images, not processed")
	parser.add_argument("-R","--redo",action="store_true",
	                    help="reprocess and overwrite existing files")
	parser.add_argument("--noclean",action="store_true",
	                    help="don't delete temporary files")
	parser.add_argument("--sdssrm",action="store_true",
	                    help="use sdssrm pipeline products")
	parser.add_argument("--naocver",action="store_true",
	                    help="use NAOC pipeline products (split by ccd)")
	parser.add_argument("-p","--processes",type=int,default=1,
	                    help="number of processes")
	args = parser.parse_args()
	if args.raw:
		process_raw_images(args.inputFiles,
		                   outputDir=args.outputdir,
		                   overwrite=args.redo,
		                   cleanup=not args.noclean)
	elif args.sdssrm:
		process_idm_images(args.inputFiles,
		                   outputDir=args.outputdir,
		                   overwrite=args.redo)
	else:
		process_images(args.inputFiles,
		               outputDir=args.outputdir,
		               overwrite=args.redo,
		               nproc=args.processes,
		               naocver=args.naocver)

