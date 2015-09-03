#!/usr/bin/env python

import os
import re
from copy import copy
from collections import OrderedDict
import numpy as np
from scipy.stats.mstats import mode
from scipy.interpolate import LSQUnivariateSpline
import fitsio
from astropy.stats import sigma_clip
from astropy.modeling import models,fitting

# the order of the amplifiers in the FITS extensions, i.e., HDU1=amp#4
ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]

# a 90Prime FITS MEF file has 16 extensions, one for each amplifier
# iterate over them in the order given above
bok90mef_extensions = ['IM%d' % a for a in ampOrder]

bokCenterAmps = ['IM4','IM7','IM10','IM13']

'''
nominal_gain =  np.array(
  [ 1.3, 1.3, 1.3, 1.3, 
    1.5, 1.5, 1.3, 1.5, 
    1.4, 1.4, 1.4, 1.3, 
    1.4, 1.3, 1.4, 1.4
   ] )
'''

nominal_gain = np.array(
  [ 1.24556017,  1.29317832,  1.31759822,  1.28293753,  
    1.44988859, 1.52633166,  1.42589855,  1.51268101,  
    1.33969975,  1.39347458, 1.3766073 ,  1.39406121,  
#    1.42733335,  1.38764536,  1.79094434, 1.45403028
    1.42733335,  1.38764536,  1.40, 1.45403028
  ] )


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

def stats_region(statreg):
	if type(statreg) is tuple:
		return statreg
	elif statreg == 'amp_central_quadrant':
		return (512,-512,512,-512)
	elif statreg == 'amp_corner_ccdcenter_small':
		return (-512,-50,-512,-50)
	elif statreg == 'amp_corner_ccdcenter':
		return (-1024,-1,-1024,-1)
	elif statreg == 'centeramp_corner_fovcenter':
		# for the 4 central amps, this is the corner towards the field center
		return (50,1024,50,1024)
	else:
		raise ValueError

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

def bok_rebin(im,nbin):
	s = np.array(im.shape) / nbin
	return im.reshape(s[0],nbin,s[1],nbin).swapaxes(1,2).reshape(s[0],s[1],-1)

def bok_getxy(hdr,coordsys='image'):
	y,x = np.indices((hdr['NAXIS2'],hdr['NAXIS1']))
	# FITS coordinates are 1-indexed (correct?)
	#x += 1
	#y += 1
	if coordsys == 'image':
		pass
	elif coordsys == 'physical':
		x = hdr['LTM1_1']*(x - hdr['LTV1'])
		y = hdr['LTM2_2']*(y - hdr['LTV2'])
	elif coordsys == 'sky':
		# hacky assumption of orthogonal coordinates but true at this stage
		dx = hdr['CD1_1'] + hdr['CD2_1']
		dy = hdr['CD1_2'] + hdr['CD2_2']
		x = np.sign(dx)*(x - hdr['CRPIX1'])
		y = np.sign(dy)*(y - hdr['CRPIX2'])
	else:
		raise ValueError
	return x,y

def bok_fov_rebin(fits,nbin,coordsys='sky'):
	rv = {}
	if type(fits) is str:
		fits = fitsio.FITS(fits)
	for hdu in fits[1:]:
		im = hdu.read()
		hdr = hdu.read_header()
		x,y = bok_getxy(hdr,coordsys)
		if nbin > 1:
			im = bok_rebin(im,nbin)
			x = x[nbin//2::nbin,nbin//2::nbin]
			y = y[nbin//2::nbin,nbin//2::nbin]
		extn = hdu.get_extname()
		rv[extn] = {'x':x,'y':y,'im':im}
	return rv

def bok_polyfit(fits,nbin,order):
	binnedIm = bok_fov_rebin(fits,nbin,'sky')
	# XXX need bad pixel and object masks here
	# collect the CCD mosaic into a single image
	X,Y,fovIm = [],[],[]
	for ccd in ['CCD%d'%i for i in range(1,5)]:
		print 'getting sky for ',ccd
		im = fits[ccd].read()
		if True:
			# hacky way of masking objects for now
			sky = sigma_clip(im[100:1500,100:1500],iters=5,sig=2.5)
			skym,skys = np.ma.median(sky),sky.std()
			mask = (im > skym+5*skys) | (im < skym-5*skys)
			print skym,skys,mask.sum()
			mask = bok_rebin(mask,nbin)
		mim = np.ma.masked_array(binnedIm[ccd]['im'],mask)
		im = mim.mean(axis=-1)
		X.append(binnedIm[ccd]['x'].flatten())
		Y.append(binnedIm[ccd]['y'].flatten())
		fovIm.append(im.flatten())
	X = np.concatenate(X)
	Y = np.concatenate(Y)
	fovIm = np.concatenate(fovIm)
	print X.shape,Y.shape,fovIm.shape
	# fit a polynomial to the binned mosaic image
	poly_model = models.Polynomial2D(degree=order)
	fitfun = fitting.LevMarLSQFitter()
	p = fitfun(poly_model,X,Y,fovIm)
	# return the model images for each CCD at native resolution
	rv = {}
	for ccd in ['CCD%d'%i for i in range(1,5)]:
		x,y = bok_getxy(fits[ccd].read_header(),'sky')
		rv[ccd] = p(x,y)
	rv['skymodel'] = p
	return rv

def bok_make_image(fits,pngfn,nbin=1,coordsys='sky',**kwargs):
	import matplotlib.pyplot as plt
	#kwargs.setdefault('cmap',plt.cm.hot_r)
	kwargs.setdefault('interpolation','nearest')
	w = 0.4575
	h = 0.465
	fov = bok_fov_rebin(fits,nbin,coordsys)
	fig = plt.figure(figsize=(6,6))
	for n,ccd in enumerate(['CCD2','CCD4','CCD1','CCD3']):
		if n == 0:
			im = fov[ccd]['im']
			i1,i2 = 100//nbin,1500//nbin
			background = sigma_clip(im[i1:i2,i1:i2],iters=3,sig=2.2)
			m,s = background.mean(),background.std()
		im = fov[ccd]['im'].mean(axis=-1)
		x = fov[ccd]['x']
		y = fov[ccd]['y']
		i = n % 2
		j = n // 2
		pos = [ 0.025 + i*w + i*0.04, (j+1)*0.025 + j*h, w, h ]
		ax = fig.add_axes(pos)
		ax.imshow(im,origin='lower',
		          extent=[x[0,0],x[0,-1],y[0,0],y[-1,0]],
		          vmin=m-2*s,vmax=m+8*s,**kwargs)
		if coordsys=='sky':
			ax.set_xlim(x.max(),x.min())
		else:
			ax.set_xlim(x.min(),x.max())
		ax.set_ylim(y.min(),y.max())

def stack_image_cube(imCube,**kwargs):
	reject = kwargs.get('reject','sigma_clip')
	method = kwargs.get('method','mean')
	scale = kwargs.get('scale')
	weights = kwargs.get('weights')
	withVariance = kwargs.get('with_variance',False)
	retScales = kwargs.get('ret_scales',False)
	x1,x2,y1,y2 = stats_region(kwargs.get('stats_region',
	                                      'amp_central_quadrant'))
	clipargs = {'iters':kwargs.get('clip_iters',2),
	            'sig':kwargs.get('clip_sig',2.5),
	            'cenfunc':np.ma.mean}
	# scale images
	if scale is not None:
		if type(scale) is np.ndarray:
			scales = scale
		elif scale.startswith('normalize'):
			imScales = imCube[y1:y2,x1:x2]/imCube[y1:y2,x1:x2,[0]]
			imScales = imScales.reshape(-1,imCube.shape[-1])
			scales = sigma_clip(imScales,cenfunc=np.ma.mean,axis=0)
			if scale.endswith('_mean'):
				scales = scales.mean(axis=0)
			else:
				# default is the mode
				scales,_ = mode(scales,axis=0)
			scales /= scales.max()
			scales **= -1
		else:
			scales = scale(imCube)
		imcube = imCube * scales
	else:
		imcube = imCube
	# mask rejected pixels
	if reject == 'sigma_clip':
		imcube = sigma_clip(imcube,axis=-1,**clipargs)
	elif reject == 'minmax':
		imcube = np.ma.masked_array(imcube)
		imcube[:,:,imcube.argmax(axis=-1)] = np.ma.masked
		imcube[:,:,imcube.argmin(axis=-1)] = np.ma.masked
	elif reject is not None:
		raise ValueError
	# do the stacking
	if method == 'mean':
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
	reject = kwargs.get('reject','sigma_clip')
	method = kwargs.get('method','mean')
	mask_at = kwargs.get('mask_at',[0,1,2,-1])
	along = kwargs.get('along','columns')
	clipargs = {'iters':kwargs.get('clip_iters',2),
	            'sig':kwargs.get('clip_sig',2.5),
	            'cenfunc':np.ma.mean}
	spline_interval = kwargs.get('spline_interval',20)
	if along == 'rows':
		# make it look like a column overscan for simplicity
		overscan = overscan.transpose()
	npix = overscan.shape[0]
	#
	overscan = np.ma.masked_array(overscan)
	overscan[:,mask_at] = np.ma.masked
	#
	if reject == 'sigma_clip':
		overscan = sigma_clip(overscan,axis=1,**clipargs)
	elif reject == 'minmax':
		overscan[:,overscan.argmax(axis=1)] = np.ma.masked
		overscan[:,overscan.argmin(axis=1)] = np.ma.masked
	#
	if method == 'mean':
		oscan_fit = overscan.mean(axis=1)
	elif method == 'mean_value':
		oscan_fit = np.repeat(overscan.mean(),npix)
	elif method == 'median_value':
		oscan_fit = np.repeat(np.ma.median(overscan),npix)
	elif method == 'cubic_spline':
		knots = np.concatenate([np.arange(1,npix,spline_interval),[npix,]])
		mean_fit = overscan.mean(axis=1)
		x = np.arange(npix)
		spl_fit = LSQUnivariateSpline(x,mean_fit,t=knots)
		oscan_fit = spl_fit(x)
	else:
		raise ValueError
	return oscan_fit

class OverscanCollection(object):
	def __init__(self,oscanImgFile,along='columns'):
		self.along = along
		self.imgFile = oscanImgFile
		if along=='columns':
			self.arr_stack = np.hstack
		else:
			self.arr_stack = np.vstack
		self.tmpfn1 = oscanImgFile+'_oscantmp.npy'
		self.tmpfn2 = oscanImgFile+'_restmp.npy'
		if os.path.exists(self.tmpfn1):
			os.unlink(self.tmpfn1)
		if os.path.exists(self.tmpfn2):
			os.unlink(self.tmpfn2)
		self.tmpOscanImgFile = open(self.tmpfn1,'ab')
		self.tmpOscanResImgFile = open(self.tmpfn2,'ab')
		self.files = []
	def close(self):
		os.unlink(self.tmpfn1)
		os.unlink(self.tmpfn2)
	def append(self,oscan,oscanFit,fileName):
		np.save(self.tmpOscanImgFile,oscan)
		if self.along=='columns':
			resim = (oscan - oscanFit[:,np.newaxis]).astype(np.float32)
		else:
			resim = (oscan - oscanFit[np.newaxis,:]).astype(np.float32)
		np.save(self.tmpOscanResImgFile,resim.filled(-999))
		self.files.append(os.path.basename(fileName))
	def write_image(self):
		self.tmpOscanImgFile.close() # could just reset file pointer?
		self.tmpOscanResImgFile.close()
		nfiles = len(self.files)
		if nfiles==0:
			return
		hdr = OrderedDict()
		hdr['NOVSCAN'] = nfiles
		for n,f in enumerate(self.files,start=1):
			hdr['OVSCN%03d'%n] = f
		f1 = open(self.tmpfn1,'rb')
		oscanImg = self.arr_stack([np.load(f1) for i in range(nfiles)])
		f1.close()
		f2 = open(self.tmpfn2,'rb')
		oscanResImg = self.arr_stack([np.load(f2) for i in range(nfiles)])
		f2.close()
		if os.path.exists(self.imgFile+'.fits'):
			os.unlink(self.imgFile+'.fits')
		oscanFits = fitsio.FITS(self.imgFile+'.fits','rw')
		oscanFits.write(oscanImg,header=hdr)
		oscanFits.write(oscanResImg,clobber=False)
		oscanFits.close()
	def n_images(self):
		return len(self.files)

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
	                              os.path.join(outputDir,'oscan_cols'))
	oscanRowsImgFile = kwargs.get('oscan_rows_file',
	                              os.path.join(outputDir,'oscan_rows'))
	if write_overscan_image:
		oscanColCollection = { extn:OverscanCollection(oscanColsImgFile+
		                                               '_'+extn)
		                          for extn in extensions }
		oscanRowCollection = { extn:OverscanCollection(oscanRowsImgFile+
		                                               '_'+extn,along='rows')
		                          for extn in extensions }
	for f in fileList:
		print '::: ',f
		fits = fitsio.FITS(f)
		outFits = fitsio.FITS(outputFileMap(f),'rw')
		hdr = fits[0].read_header()
		outFits.write(None,header=hdr)
		for extn in extensions:
			data,oscan_cols,oscan_rows = extract_overscan(fits[extn])
			colbias = fit_overscan(oscan_cols,**kwargs)
			data[:] -= colbias[:,np.newaxis]
			if oscan_rows is not None:
				# first fit and then subtract the overscan columns at the
				# end of the strip of overscan rows
				# XXX hardcoded geometry
				_colbias = fit_overscan(oscan_rows[:,-20:],**kwargs)
				oscan_rows = oscan_rows[:,:-20] - _colbias[:,np.newaxis]
				# now fit and subtract the overscan rows
				rowbias = fit_overscan(oscan_rows,along='rows',**kwargs)
				data[:] -= rowbias[np.newaxis,:data.shape[1]]
			# write the output file
			hdr = fits[extn].read_header()
			hdr['OSCANSUB'] = 'method=%s' % kwargs.get('method','default')
			# something changed about what median returns...
			try:
				hdr['OSCANMED'] = float(np.ma.median(colbias).filled(-999))
			except:
				hdr['OSCANMED'] = float(np.ma.median(colbias))
			outFits.write(data,extname=extn,header=hdr)
			# save the oscans to images
			if write_overscan_image:
				oscanColCollection[extn].append(oscan_cols,colbias,f)
				if oscan_rows is not None:
					oscanRowCollection[extn].append(oscan_rows,rowbias,f)
		fits.close()
		outFits.close()
	if write_overscan_image:
		for extn in extensions:
			oscanColCollection[extn].write_image()
			oscanColCollection[extn].close()
			if oscanRowCollection[extn].n_images() > 0:
				oscanRowCollection[extn].write_image()
			oscanRowCollection[extn].close()

def stack_bias_frames(fileList,**kwargs):
	extensions = kwargs.get('extensions',bok90mef_extensions)
	outputDir = kwargs.get('output_dir','./')
	outputFile = kwargs.get('output_file','bias.fits')
	withVariance = kwargs.get('with_variance',False)
	varOutputFile = kwargs.get('var_output_file','biasvar.fits')
	kwargs.setdefault('method','clipped_mean')
	kwargs.setdefault('clip_iters',1)
	check_dropped_rows = kwargs.get('check_dropped_rows',False)
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
	#doIllumCorr = kwargs.get('illum_corr',True)
	withVariance = kwargs.get('with_variance',False)
	varOutputFile = kwargs.get('var_output_file','biasvar.fits')
	retainCounts = kwargs.get('retain_counts',False)
	x1,x2,y1,y2 = stats_region(kwargs.get('stats_region',
	                                      'amp_corner_ccdcenter_small'))
	                                      #'amp_central_quadrant'))
	_kwargs = copy(kwargs)
	_kwargs.setdefault('scale','normalize')
	_kwargs.setdefault('ret_scales',True)
	if biasFile is None:
		biasFits = None
	else:
		biasFits = fitsio.FITS(biasFile)
	fits = fitsio.FITS(os.path.join(outputDir,outputFile),'rw')
	hdr = _write_stack_header_cards(fileList,'FLAT')
	fits.write(None,header=hdr)
	if withVariance:
		varFits = fitsio.FITS(os.path.join(outputDir,varOutputFile),'rw')
		varFits.write(None,header=hdr)
	for extn in extensions:
		print '::: %s extn %s' % (outputFile,extn)
		flatCube = build_cube(fileList,extn)
		if biasFits is not None:
			flatCube -= biasFits[extn].read()[:,:,np.newaxis]
		stack,extras = stack_image_cube(flatCube,**_kwargs)
		if retainCounts:
			flatNorm = 1.0  # do not normalize to unity
		else:
			flatNorm = mode(stack[y1:y2,x1:x2],axis=None)[0][0]
		stack /= flatNorm
		stack = stack.filled(1.0)
		hdr = fitsio.read_header(fileList[0],extn)
		scales = extras[0].squeeze().filled(np.nan)
		hdr['FLATNORM'] = flatNorm
		for _i,_scl in enumerate(scales,start=1):
			hdr['FLTSCL%02d'%_i] = _scl
		fits.write(stack,extname=extn,header=hdr)
		if withVariance:
			varFits.write(extras[0],extname=extn,header=hdr)
	fits.close()
	if biasFits is not None:
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
	hdr = _write_stack_header_cards(fileList,'SKYFLAT')
	fits.write(None,header=hdr)
	for extn in extensions:
		flatCube = build_cube(fileList,extn,masks=objectMasks)
		stack,extras = stack_image_cube(flatCube,**_kwargs)
		flatNorm = mode(stack[y1:y2,x1:x2],axis=None)[0][0]
		stack /= flatNorm
		# XXX smooth it
		stack = stack.filled(1.0)
		hdr = fitsio.read_header(fileList[0],extn)
		hdr['FLATNORM'] = flatNorm
		fits.write(stack,extname=extn,header=hdr)
	fits.close()



###############################################################################
#                                                                             #
#                              PROCESS ROUND 1                                #
#                   bias subtraction & flat field                             #
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
		print f
		if outputFileMap is None:
			# modify the file in-place
			outFits = fitsio.FITS(f,'rw')
			inFits = outFits
			hdr0 = inFits[0].read_header()
		else:
			inFits = fitsio.FITS(f)
			hdr0 = inFits[0].read_header()
			outFits = fitsio.FITS(outputFileMap(f),'rw')
			outFits.write(None,header=hdr0)
		if biasSubMap is not None:
			biasSubFits = fitsio.FITS(biasSubMap(f),'rw')
			biasSubFits.write(None,header=hdr0)
		if flatDivMap is not None:
			flatDivFits = fitsio.FITS(flatDivMap(f),'rw')
			flatDivFits.write(None,header=hdr0)
		for extn in extensions:
			data = inFits[extn].read()
			hdr = inFits[extn].read_header()
			# BIAS SUBTRACTION
			data -= biasFits[extn][:,:]
			if biasSubMap is not None:
				biasSubFits.write(data,extname=extn,header=hdr)
			# FIRST ORDER FLAT FIELD CORRECTION
			data /= flatFits[extn][:,:]
			if flatDivMap is not None:
				flatDivFits.write(data,extname=extn,header=hdr)
			hdr['BIASFILE'] = biasFile
			hdr['FLATFILE'] = flatFile
			outFits.write(data,extname=extn,header=hdr)
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
#                balances the amplifiers with a gain correction               #
#                                                                             #
###############################################################################

def multiply_gain(inFits,extGroup,hdr,skyGainCor,inputGain,
                  ampCorStatReg,ccdCorStatReg,clipArgs,skyIn,refAmp):
	# the mode doesn't seem to be well-behaved here (why?), 
	#sky_est = lambda x: mode(x,axis=None)[0][0]
	# mean seems robust
	#sky_est = np.ma.median
	sky_est = np.ma.mean
	# the stats region used to balance amps within a CCD using sky values
	xa1,xa2,ya1,ya2 = ampCorStatReg
	# the stats region used to balance CCDs across the field using sky values
	xc1,xc2,yc1,yc2 = ccdCorStatReg
	# load the per-amp images
	ampIms = [ inFits[ext].read() for ext in extGroup ]
	# start with the input gain values (from header keywords or input by user)
	gain = np.array([ inputGain[ext] for ext in extGroup ])
	# use the sky counts to balance the gains
	if skyGainCor:
		# first balance across amplifers
		rawSky = np.array([ sky_est(sigma_clip(im[ya1:ya2,xa1:xa2],
		                                       **clipArgs))
		                           for im in ampIms ])
		skyCounts = rawSky * gain
		refAmpIndex = np.where(extGroup == refAmp)[0][0]
		gain2 = skyCounts[refAmpIndex] / skyCounts
		# then balance across CCDs
		centerAmp = (set(extGroup) & set(bokCenterAmps)).pop()
		ci = np.where(extGroup == centerAmp)[0][0]
		skyCounts = sky_est(sigma_clip(ampIms[ci][yc1:yc2,xc1:xc2],**clipArgs))
		__rawSky2 = skyCounts
		skyCounts *= gain[ci] * gain2[ci]
		if skyIn is None:
			gain3 = 1.0
		else:
			gain3 = skyIn / skyCounts
		for i,ext in enumerate(extGroup):
			chNum = int(ext.replace('IM',''))
			hdr['SKYC%02d%s1'%(chNum,'ABCD'[i])] = rawSky[i]
			hdr['GAIN%02d%s1'%(chNum,'ABCD'[i])] = gain[i]
			hdr['GAIN%02d%s2'%(chNum,'ABCD'[i])] = gain2[i]
		hdr['CCDGAIN3'] = gain3
		gain *= gain2 * gain3
	else:
		skyCounts = None
		# store the (default) gain values used
		for i,ext in enumerate(extGroup):
			chNum = int(ext.replace('IM',''))
			hdr['GAIN%02d%s1'%(chNum,'ABCD'[i])] = gain[i]
	ampIms = [ im*g for im,g in zip(ampIms,gain) ]
	return ampIms,skyCounts

def _orient_mosaic(hdr,ims,ccdNum,origin):
	im1,im2,im3,im4 = ims
	# orient the channel images N through E and stack into CCD image
	outIm = np.vstack([ np.hstack([ np.flipud(np.rot90(im2)),
	                                np.rot90(im4,3) ]),
	                    np.hstack([ np.rot90(im1),
	                                np.fliplr(np.rot90(im3)) ]) ])
	if origin == 'lower left':
		pass
	elif origin == 'center':
		if ccdNum == 1:
			outIm = np.fliplr(outIm)
		elif ccdNum == 2:
			outIm = np.rot90(outIm,2)
		elif ccdNum == 3:
			pass
		elif ccdNum == 4:
			outIm = np.flipud(outIm)
	ny,nx = outIm.shape
	det_i = (ccdNum-1) // 2
	det_j = ccdNum % 2
	hdr['DATASEC'] = '[1:%d,1:%d]' % (nx,ny)
	hdr['DETSEC'] = '[%d:%d,%d:%d]' % (nx*det_i+1,nx*(det_i+1),
	                                   ny*det_j+1,ny*(det_j+1))
	plateScale = np.max(np.abs([hdr['CD1_1'],hdr['CD1_2']]))
	# --> these two disagree in physical coords by 1 pixel for det_i==1
	if origin == 'center':
		# works for WCS but not IRAF for some reason
		hdr['CD1_1'] = 0.0
		hdr['CD2_2'] = 0.0
		signx = [-1,+1][det_i]
		signy = [-1,+1][det_j]
		hdr['CD2_1'] = -signx*plateScale
		hdr['CD1_2'] = signy*plateScale
		hdr['CRPIX1'] = -182.01
		hdr['CRPIX2'] = -59.04
		hdr['LTM1_1'] = float(signx)
		hdr['LTM2_2'] = float(signy)
		hdr['LTV1'] = [4096.0,-4097.0][det_i]
		hdr['LTV2'] = [4033.0,-4032.0][det_j]
	elif origin == 'lower left':
		hdr['LTM1_1'] = 1.0
		hdr['LTM2_2'] = 1.0
		hdr['LTV1'] = 0 if det_i == 0 else -nx
		hdr['LTV2'] = 0 if det_j == 0 else -ny
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
	return outIm,hdr

def combine_ccds(fileList,**kwargs):
	outputFileMap = kwargs.get('output_file_map')
	tmpFileName = 'tmp.fits'
	# do the extensions in numerical order, instead of HDU list order
	extns = np.array(['IM%d' % ampNum for ampNum in range(1,17)])
	#
	inputGain = kwargs.get('input_gain')
	skyGainCor = kwargs.get('sky_gain_correct',True)
	ampCorStatReg = stats_region(kwargs.get('stats_region',
	                                      'amp_corner_ccdcenter'))
	# not a keyword (?)
	ccdCorStatReg = stats_region('centeramp_corner_fovcenter')
	clipArgs = {'iters':kwargs.get('clip_iters',3),
	            'sig':kwargs.get('clip_sig',2.5),
	            'cenfunc':np.ma.mean}
	#origin = kwargs.get('origin','lower left')
	origin = kwargs.get('origin','center')
	# 2,8 are fairly stable and not in corner (4 is worst on CCD1, 7 on CCD2)
	# 11 is on CCD3 but not affected by bias ramp
	# 16 is by far the least affected by A/D errors on CCD4
	refAmps = ['IM2','IM8','IM11','IM16']
	#refAmps = ['IM4','IM8','IM9','IM13']
	if inputGain is None:
		inputGain = { 'IM%d'%ampNum:g 
		                  for ampNum,g in zip(ampOrder,nominal_gain)}
	#
	for f in fileList:
		print 'combine: ',f
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
		refSkyCounts = None
		for ccdNum,extGroup in enumerate(np.hsplit(extns,4),start=1):
			hdr = inFits[bokCenterAmps[ccdNum-1]].read_header()
			# load the individual channels and balance them with a gain
			# correction (either default values or using sky counts)
			(im1,im2,im3,im4),skyCounts = \
			        multiply_gain(inFits,extGroup,hdr,skyGainCor,
			                      inputGain,ampCorStatReg,ccdCorStatReg,
			                      clipArgs,refSkyCounts,refAmps[ccdNum-1])
			if ccdNum == 1:
				refSkyCounts = skyCounts
			# orient the channel images into a mosaic of CCDs and
			# modify WCS & mosaic keywords
			outIm,hdr = _orient_mosaic(hdr,(im1,im2,im3,im4),ccdNum,origin)
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

def subtract_sky(fileList,**kwargs):
	outputFileMap = kwargs.get('output_file_map')
	for f in fileList:
		print f
		if outputFileMap is None:
			# modify the file in-place
			outFits = fitsio.FITS(f,'rw')
			inFits = outFits
			hdr0 = inFits[0].read_header()
		else:
			inFits = fitsio.FITS(f)
			hdr0 = inFits[0].read_header()
			outFits = fitsio.FITS(outputFileMap(f),'rw')
			outFits.write(None,header=hdr0)
		#
		skyFit = bok_polyfit(inFits,64,1)
		# subtract the sky level at the origin so as to only remove 
		# the gradient
		sky0 = skyFit['skymodel'](0,0)
		for ccd in ['CCD%d'%i for i in range(1,5)]:
			hdr = inFits[ccd].read_header()
			data = inFits[ccd][:,:] - (skyFit[ccd] - sky0)
			hdr['SKYVAL'] = sky0
			outFits.write(data,extname=ccd,header=hdr)
		if outFits != inFits:
			outFits.close()
		inFits.close()

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



