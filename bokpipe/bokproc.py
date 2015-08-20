#!/usr/bin/env python

import os
import re
from copy import copy
import numpy as np
import fitsio
from astropy.stats import sigma_clip
from scipy.stats.mstats import mode

# the order of the amplifiers in the FITS extensions, i.e., HDU1=amp#4
ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]

# a 90Prime FITS MEF file has 16 extensions, one for each amplifier
# iterate over them in the order given above
bok90mef_extensions = ['IM%d' % a for a in ampOrder]


###############################################################################
#                                                                             #
#                          GENERAL UTILITIES                                  #
#                                                                             #
###############################################################################

class OutputExistsError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class FileNameMap(object):
	def __init__(self,newDir=None,newSuffix=None,strip_gz=True):
		self.newDir = newDir
		self.newSuffix = newSuffix
		self.strip_gz = strip_gz
	def __call__(self,fileName):
		if self.newDir is None:
			newDir = os.path.dirname(fileName)
		else:
			newDir = self.newDir
		fn = os.path.basename(fileName)
		if self.strip_gz and fn.endswith('.gz'):
			fn = fn[:-3]
		if self.newSuffix is not None:
			fn = fn.replace('.fits',self.newSuffix+'.fits')
		return os.path.join(newDir,fn)

def _convertfitsreg(regstr):
	regpattern = r'\[(\d+):(\d+),(\d+):(\d+)\]'
	rv =  [ int(d) for d in  re.match(regpattern,regstr).groups() ]
	# FITS region indices are 1-indexed
	rv[0] -= 1
	rv[2] -= 1
	return rv

def build_cube(fileList,extn,masks=None):
	cube = np.dstack( [ fitsio.read(f,extn) for f in fileList ] )
	if masks is not None:
		if isinstance(masks,FileNameMap):
			mask = np.dstack([ fitsio.read(masks(f),extn) for f in fileList ])
		else:
			mask = np.dstack([ fitsio.read(f,extn) for f in masks ])
	else:
		mask = None
	cube = np.ma.masked_array(cube,mask)
	return cube

def stack_image_cube(imCube,**kwargs):
	method = kwargs.get('method','clipped_mean')
	scale = kwargs.get('scale')
	weights = kwargs.get('weights')
	withVariance = kwargs.get('with_variance',False)
	retScales = kwargs.get('ret_scales',False)
	statRegion = kwargs.get('stat_region',(512,-512,512,-512))
	clipargs = {'iters':kwargs.get('clip_iters',2),
	            'sig':kwargs.get('clip_sig',2.5)}
	x1,x2,y1,y2 = statRegion
	# scale images
	if scale is not None:
		if type(scale) is np.ndarray:
			scales = scale
		elif scale.startswith('normalize'):
			imScales = imCube[y1:y2,x1:x2]/imCube[y1:y2,x1:x2,[0]]
			imScales = imScales.reshape(-1,imCube.shape[-1])
			scales = sigma_clip(imScales,cenfunc=np.mean,axis=0)
			if scale.endswith('_mean'):
				scales = scales.mean(axis=0)
			else:
				# default is the mode
				scales,_ = mode(scales,axis=0)
			scales /= scales.max()
		else:
			scales = scale(imCube)
		imcube = imCube * scales
	else:
		imcube = imCube
	# do the stacking
	if method == 'clipped_mean':
		imcube = sigma_clip(imcube,axis=-1,**clipargs)
		stack = np.ma.average(imcube,axis=-1,weights=weights)
	elif method == 'mean':
		stack = np.ma.average(imcube,weights=weights,axis=-1)
	elif method == 'median':
		stack = np.ma.median(imcube,axis=-1)
	else:
		raise ValueError
	# why does it get upcasted to float64?
	stack = stack.astype(np.float32)
	extras = []
	if retScales:
		extras.append(scales)
	if withVariance:
		var = np.ma.var(imcube,axis=-1).filled(0).astype(np.float32)
		extras.append(var)
	return stack,extras

def _write_stack_header_cards(fileList,cardPrefix):
	hdr = fitsio.read_header(fileList[0])
	for num,f in enumerate(fileList,start=1):
		hdr['%s%03d'%(cardPrefix,num)] = os.path.basename(f)
	return hdr




###############################################################################
#                                                                             #
#                            BIAS SUBTRACTION                                 #
#                                                                             #
###############################################################################

def extract_overscan(imhdu):
	'''Given a 90Prime FITS HDU corresponding to a single amplifier, extract
	   the overscan regions and trim the image.
	   Returns data, overscan_cols, overscan_rows
	   Output is converted to floats
	'''
	data = imhdu.read()
	hdr = imhdu.read_header()
	x1,x2,y1,y2 = _convertfitsreg(hdr['BIASSEC'])
	overscan_cols = data[y1:y2,x1:x2].astype(np.float32)
	x1,x2,y1,y2 = _convertfitsreg(hdr['DATASEC'])
	if hdr['NAXIS2'] > y2+1:
		# added overscan rows are not identified in header keywords, just
		# extract any extra rows outside of DATASEC
		overscan_rows = data[y2:,:].astype(np.float32)
	else:
		overscan_rows = None
	data = data[y1:y2,x1:x2].astype(np.float32)
	return ( data,overscan_cols,overscan_rows )

def fit_overscan(overscan,**kwargs):
	method = kwargs.get('method','clipped_mean')
	mask_columns = kwargs.get('mask_columns',[0,1,-1])
	along = kwargs.get('along','columns')
	clipargs = {'iters':kwargs.get('clip_iters',2),
	            'sig':kwargs.get('clip_sig',2.5)}
	if along == 'columns':
		axis = 1
	elif along == 'rows':
		axis = 0
	else:
		raise ValueError
	#
	mask = np.zeros(overscan.shape,dtype=bool)
	if mask_columns is not None:
		mask[:,mask_columns] = True
	overscan = np.ma.masked_array(overscan,mask=mask)
	#
	if method == 'clipped_mean':
		oscan_fit = sigma_clip(overscan,axis=axis,**clipargs).mean(axis=axis)
	elif method == 'clipped_mean_value':
		oscan_fit = sigma_clip(overscan,**clipargs).mean(axis=axis)
	elif method == 'median_value':
		oscan_fit = np.repeat(np.ma.median(overscan),overscan.shape[0])
	else:
		raise ValueError
	return oscan_fit

def _write_oscan_image(f,oscanImgFile,oscan,oscanFit,along='columns'):
	arr_stack = np.hstack if along=='columns' else np.vstack
	oscanFits = fitsio.FITS(oscanImgFile,'rw')
	rv = oscanFits[0].read(header=True)
	if rv is None:
		oscanFits.write(oscanImg,oscan,clobber=True,
		                header={'NOVSCAN':1,'OVSCN001':os.path.basename(f)})
		oscanFits.write(oscanResImg,oscan-oscanFit)
	else:
		oscanImg,hdr = rv
		hdr['NOVSCAN'] += 1
		hdr['OVSCN%03d'%hdr['NOVSCAN']] = os.path.basename(f)
		oscanFits[0].write(arr_stack([oscanImg,oscan]),header=hdr)
		oscanResImg,hdr = oscanFits[1].read()
		oscanFits[1].write(arr_stack([oscanResImg,oscan-oscanFit]))
	oscanFits.close()

def _imsub(f1,f2,outf,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	fits1 = fitsio.FITS(f1)
	fits2 = fitsio.FITS(f2)
	outFits = fitsio.FITS(outf,'rw')
	hdr = fits1[0].read_header()
	outFits.write(None,header=hdr)
	for extn in extensions:
		data = fits1[extn][:,:] - fits2[extn][:,:]
		hdr = fits1[extn].read_header()
		outFits.write(data,extname=extn,header=hdr)
	fits1.close()
	fits2.close()
	outFits.close()

def subtract_overscan(fileList,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	# XXX needs to have a default
	outputFileMap = kwargs.get('output_file_map')
	outputDir = './'
	write_overscan_image = kwargs.get('write_overscan_image',False)
	oscanColsImgFile = kwargs.get('oscan_cols_file',
	                              os.path.join(outputDir,'oscan_cols.fits'))
	oscanRowsImgFile = kwargs.get('oscan_rows_file',
	                              os.path.join(outputDir,'oscan_rows.fits'))
	for f in fileList:
		fits = fitsio.FITS(f)
		outFits = fitsio.FITS(outputFileMap(f),'rw')
		hdr = fits[0].read_header()
		outFits.write(None,header=hdr)
		for extn in extensions:
			data,oscan_cols,oscan_rows = extract_overscan(fits[extn])
			colbias = fit_overscan(oscan_cols,**kwargs)
			data[:] -= colbias[:,np.newaxis]
			if oscan_rows is not None:
				rowbias = fit_overscan(oscan_rows,along='rows',**kwargs)
				# have to trim overscan edge
				data[:] -= rowbias[np.newaxis,:data.shape[1]]
			# write the output file
			hdr = fits[extn].read_header()
			hdr['OSCANSUB'] = 'method=%s' % kwargs.get('method','default')
			hdr['OSCANMED'] = np.median(colbias)
			outFits.write(data,extname=extn,header=hdr)
			# save the oscans to images
			if write_overscan_image:
				_write_oscan_image(f,oscanColsImgFile,oscan_cols,colbias)
				if oscan_rows is not None:
					_write_oscan_image(f,oscanRowsImgFile,oscan_rows,rowbias,
					                   along='rows')
		fits.close()
		outFits.close()

def stack_bias_frames(fileList,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	outputDir = kwargs.get('output_dir','./')
	outputFile = kwargs.get('output_file','bias.fits')
	withVariance = kwargs.get('with_variance',False)
	varOutputFile = kwargs.get('var_output_file','biasvar.fits')
	kwargs.setdefault('method','mean')
	check_dropped_rows = kwargs.get('check_dropped_rows',True)
	fits = fitsio.FITS(os.path.join(outputDir,outputFile),'rw')
	hdr = _write_stack_header_cards(fileList,'BIAS')
	fits.write(None,header=hdr)
	if withVariance:
		varFits = fitsio.FITS(os.path.join(outputDir,varOutputFile),'rw')
		varFits.write(None,header=hdr)
	for extn in extensions:
		print '::: %s extn %s' % (outputFile,extn)
		biasCube = build_cube(fileList,extn)
		if check_dropped_rows:
			# slice out columns from the center of the array
			centerCol = biasCube.shape[1]/2
			centerSlices = biasCube[:,centerCol-10:centerCol+11]
			for k in range(biasCube.shape[-1]):
				# obtain a median column vector from the central region
				medianCol = np.median(centerSlices[:,:,k],axis=1)
				# get the pixel statistics from the middle rows
				s = sigma_clip(centerSlices[200:-200])
				# look for aberrant rows at the bottom of the array
				rej = np.abs(medianCol[:200]-s.mean())/s.std() > 3.0
				# find the first row not rejected
				row1 = np.where(~rej)[0][0]
				# if more than 10 rows are bad, must be a bad image
				if row1 > 10:
					# ... and pad the mask a bit
					row1 += 3
				# now mask the affected rows
				biasCube[:row1,:,k] = np.ma.masked
		stack,extras = stack_image_cube(biasCube,**kwargs)
		# handles the dropped rows if present
		if stack.mask.any():
			rowMedian = np.ma.median(stack,axis=0)
			rowFill = np.tile(rowMedian,(stack.shape[0],1))
			stack[stack.mask] = rowFill[stack.mask]
			stack.mask[:] = False
		stack = stack.filled()
		hdr = fitsio.read_header(fileList[0],extn)
		fits.write(stack,extname=extn,header=hdr)
		if withVariance:
			varFits.write(extras[0],extname=extn,header=hdr)
	fits.close()
	if withVariance:
		varFits.close()




###############################################################################
#                                                                             #
#                                FLAT FIELDS                                  #
#                                                                             #
###############################################################################

def stack_flat_frames(fileList,biasFile,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	outputDir = kwargs.get('output_dir','./')
	outputFile = kwargs.get('output_file','flat.fits')
	withVariance = kwargs.get('with_variance',False)
	varOutputFile = kwargs.get('var_output_file','biasvar.fits')
	statRegion = kwargs.get('stat_region',(512,-512,512,-512))
	x1,x2,y1,y2 = statRegion
	_kwargs = copy(kwargs)
	_kwargs.setdefault('scale','normalize')
	biasFits = fitsio.FITS(biasFile)
	fits = fitsio.FITS(os.path.join(outputDir,outputFile),'rw')
	hdr = _write_stack_header_cards(fileList,'BIAS')
	fits.write(None,header=hdr)
	if withVariance:
		varFits = fitsio.FITS(os.path.join(outputDir,varOutputFile),'rw')
		varFits.write(None,header=hdr)
	for extn in extensions:
		print '::: %s extn %s' % (outputFile,extn)
		flatCube = build_cube(fileList,extn)
		flatCube -= biasFits[extn].read()[:,:,np.newaxis]
		stack,extras = stack_image_cube(flatCube,**_kwargs)
		stack /= mode(stack[y1:y2,x1:x2],axis=None)[0]
		stack = stack.filled(1.0)
		hdr = fitsio.read_header(fileList[0],extn)
		fits.write(stack,extname=extn,header=hdr)
		if withVariance:
			varFits.write(extras[0],extname=extn,header=hdr)
	fits.close()
	biasFits.close()
	if withVariance:
		varFits.close()

# XXX num running images
def make_supersky_flats(fileList,objectMasks,**kwargs):
	#extensions = kwargs.get('extensions',bok90mef_extensions)
	extensions = ['CCD%d' % i for i in range(1,5)]
	outputDir = kwargs.get('output_dir','./')
	outputFile = kwargs.get('output_file','superskyflat.fits')
	_kwargs = copy(kwargs)
	_kwargs.setdefault('scale','normalize')
	fits = fitsio.FITS(os.path.join(outputDir,outputFile),'rw')
	for extn in extensions:
		flatCube = build_cube(fileList,extn,masks=objectMasks)
		stack,extras = stack_image_cube(flatCube,**_kwargs)
		stack /= mode(stack[y1:y2,x1:x2])
		# XXX smooth it
		stack = stack.filled(1.0)
		hdr = _write_stack_header_cards(fileList,'FLAT')
		fits.write(stack,extname=extn,header=hdr)
	fits.close()



###############################################################################
#                                                                             #
#                              PROCESS ROUND 1                                #
#                   bias subtraction, flat field, gain multiply               #
#                                                                             #
###############################################################################

def process_round1(fileList,biasFile,flatFile,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	outputFileMap = kwargs.get('output_file_map')
	biasSubMap = kwargs.get('bias_sub_map')
	flatDivMap = kwargs.get('flat_div_map')
	biasFits = fitsio.FITS(biasFile)
	flatFits = fitsio.FITS(flatFile)
	for f in fileList:
		if outputFileMap is None:
			# modify the file in-place
			outFits = fitsio.FITS(f,'rw')
			inFits = outFits
		else:
			inFits = fitsio.FITS(f)
			outFits = fitsio.FITS(outputFileMap(f),'rw')
		if biasSubMap is not None:
			biasSubFits = fitsio.FITS(biasSubMap(f),'rw')
		if flatDivMap is not None:
			flatDivFits = fitsio.FITS(flatDivMap(f),'rw')
		meanSky = []
		for extn in extensions:
			data,hdr = inFits[extn].read(header=True)
			data -= biasFits[extn][:,:]
			if biasSubMap is not None:
				biasSubFits.write(data,extname=extn,header=hdr)
			data /= flatFits[extn][:,:]
			if flatDivMap is not None:
				flatDivFits.write(data,extname=extn,header=hdr)
			hdr['BIASFILE'] = biasFile
			hdr['FLATFILE'] = flatFile
			data *= nominal_gain[extn]
			outFits.write(data,extname=extn,header=hdr)
			sky = sigma_clip(data[y1:y2,x1:x2],iters=5,sig=2.5)
			meanSky.append(mode(sky))
		# rescale gain by sky background
		meanSky = np.array(meanSky)
		gainCorrect = meanSky[0] / meanSky
		for i,extn in enumerate(extensions):
			data,hdr = outFits[extn].read(header=True)
			g = nominal_gain[extn] * gainCorrect[i]
			data *= g
			hdr['SKYADU'] = meanSky[i]
			hdr['GAIN1'] = nominal_gain[extn]
			hdr['GAIN2'] = gainCorrect[i]
			outFits[extn].write(data,header=hdr)
		if outFits != inFits:
			outFits.close()
		inFits.close()
		if biasSubMap is not None:
			biasSubFits.close()
		if flatDivMap is not None:
			flatDivFits.close()
	biasFits.close()
	flatFits.close()



###############################################################################
#                                                                             #
#                               COMBINE CCDs                                  #
#                                                                             #
###############################################################################

def combine_ccds(fileList,**kwargs):
	outputFileMap = kwargs.get('output_file_map')
	tmpFileName = 'tmp.fits'
	# do the extensions in numerical order, instead of HDU list order
	extns = np.array(['IM%d' % ampNum for ampNum in range(1,17)])
	centerAmp = {1:'IM4',2:'IM7',3:'IM10',4:'IM13'}
	for f in fileList:
		inFits = fitsio.FITS(f)
		if outputFileMap is not None:
			outFits = fitsio.FITS(outputFileMap(f),'rw')
		else:
			# have to use a temporary file to change format
			outFits = fitsio.FITS(tmpFileName,'rw')
		hdr = inFits[0].read_header()
		hdr['DETSIZE'] = '[1:%d,1:%d]' % (8192,8064) # hardcoded
		hdr['NEXTEND'] = 4
		outFits.write(None,header=hdr)
		for ccdNum,extGroup in enumerate(np.hsplit(extns,4),start=1):
			im1,im2,im3,im4 = [inFits[ext].read() for ext in extGroup]
			hdr = inFits[centerAmp[ccdNum]].read_header()
			outIm = np.vstack([ np.hstack([ np.flipud(np.rot90(im2)),
			                                np.rot90(im4,3) ]),
			                    np.hstack([ np.rot90(im1),
			                                np.fliplr(np.rot90(im3)) ]) ])
			# modify WCS & mosaic keywords
			ny,nx = outIm.shape
			det_i = (ccdNum-1) // 2
			det_j = ccdNum % 2
			hdr['DATASEC'] = '[1:%d,1:%d]' % (nx,ny)
			hdr['DETSEC'] = '[%d:%d,%d:%d]' % (nx*det_i+1,nx*(det_i+1),
			                                   ny*det_j+1,ny*(det_j+1))
			hdr['LTM1_1'] = 1.0
			hdr['LTM2_2'] = 1.0
			hdr['LTV1'] = 0 if det_i == 0 else -nx
			hdr['LTV2'] = 0 if det_j == 0 else -ny
			plateScale = np.max(np.abs([hdr['CD1_1'],hdr['CD1_2']]))
			hdr['CD1_1'] = 0.0
			hdr['CD1_2'] = plateScale
			hdr['CD2_1'] = -plateScale
			hdr['CD2_2'] = 0.0
			crpix1,crpix2 = hdr['CRPIX1'],hdr['CRPIX2']
			if det_i==0:
				hdr['CRPIX1'] = 1 + nx - crpix2  # not really sure why +1
			else:
				hdr['CRPIX1'] = 1 + crpix2
			if det_j==0:
				hdr['CRPIX2'] = ny - crpix1
			else:
				hdr['CRPIX2'] = crpix1
			outFits.write(outIm,extname='CCD%d'%ccdNum,header=hdr)
		outFits.close()
		if outputFileMap is None:
			os.rename(tmpFileName,f)



###############################################################################
#                                                                             #
#                              PROCESS ROUND 2                                #
#   obj detection & masking, divide by supersky flat, ... (CR rejection?)     #
#                                                                             #
###############################################################################

def sextract_pass1(fileList,**kwargs):
	catalogFileNameMap = kwargs.get('catalog_name_map',FileNameMap('.cat'))
	withPsf = kwargs.get('with_psf',False)
	objMaskFileMap = kwargs.get('object_mask_map')
	for f in fileList:
		catalogFile = catalogFileNameMap(f)
		cmd = ['sex','-c','config/bok_pass1.sex',
		       '-CATALOG_NAME',catalogFile]
		if objMaskFileMap:
			cmd = cmd.extend(['-CHECKIMAGE_TYPE','OBJECT',
			                  '-CHECKIMAGE_NAME',objMaskFileMap(f)])
		cmd = cmd.append(f)
		subprocess.call(cmd)

def process_round2(fileList,superSkyFlatFile,**kwargs):
	outputFileMap = kwargs.get('output_file_map')
	flatDivMap = kwargs.get('flat_div_map')
	skyFlatFits = fitsio.FITS(superSkyFlatFile)
	sextract_pass1(fileList,**kwargs)
	make_supersky_flats(fileList,**kwargs)
	#extensions = kwargs.get('extensions',bok90mef_extensions)
	extensions = ['CCD%d' % i for i in range(1,5)]
	for f in fileList:
		if outputFileMap is None:
			# modify the file in-place
			outFits = fitsio.FITS(f,'rw')
			inFits = outFits
		else:
			inFits = fitsio.FITS(f)
			outFits = fitsio.FITS(outputFileMap(f),'rw')
		if flatDivMap is not None:
			flatDivFits = fitsio.FITS(flatDivMap(f),'rw')
		for extn in extensions:
			data,hdr = inFits[extn].read(header=True)
			data /= skyFlatFits[extn][:,:]
			if flatDivMap is not None:
				flatDivFits.write(data,extname=extn,header=hdr)
			# XXX now gain-correct using the sky
			hdr['SKYFLATF'] = superSkyFlatFile
			outFits.write(data,extname=extn,header=hdr)
		if outFits != inFits:
			outFits.close()
		inFits.close()
		if flatDivMap is not None:
			flatDivFits.close()
	skyFlatFits.close()



