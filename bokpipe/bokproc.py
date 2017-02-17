#!/usr/bin/env python

import os,sys
import re
import shutil
import tempfile
import multiprocessing
from functools import partial
import numpy as np
from scipy.interpolate import LSQBivariateSpline,RectBivariateSpline,griddata
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import interp1d,interp2d
from scipy.signal import spline_filter
from scipy.ndimage.morphology import binary_dilation,binary_closing
import scipy.ndimage.measurements as meas
from skimage import measure
from astropy.stats import sigma_clip
from astropy.modeling import models,fitting
from astropy.convolution.convolve import convolve
from astropy.convolution.kernels import Gaussian2DKernel
import fitsio

from .bokio import *
from . import bokdm
from . import bokutil
from .bokoscan import extract_overscan,overscan_subtract

# the order of the amplifiers in the FITS extensions, i.e., HDU1=amp#4
ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]
def amp_iterator():
	'''iterate over amplifiers split by CCD'''
	for i in range(0,16,4):
		yield ampOrder[i:i+4]

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
    1.42733335,  1.38764536,  1.3538, 1.45403028
  ] )

saturation_dn = 62000

# XXX
configdir = os.environ['BOKPIPE']+'/config/'

class OverExposedFlatException(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class BokImArith(bokutil.BokProcess):
	def __init__(self,op,operand,**kwargs):
		super(BokImArith,self).__init__(**kwargs)
		ops = {'+':np.add,'-':np.subtract,'*':np.multiply,'/':np.divide}
		try:
			self.op = ops[op]
		except:
			raise ValueError("operation %s not supported" % op)
		self.operand = operand
		self.operandFits = fitsio.FITS(self.operand)
	def process_hdu(self,extName,data,hdr):
		return self.op(data,self.operandFits[extName][:,:]),hdr

class BokImStat(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokImStat,self).__init__(**kwargs)
		self.statSec = bokutil.stats_region(kwargs.get('stats_region'),
		                                    kwargs.get('stats_stride'))
		self.clipArgs = kwargs.get('clip_args',{})
		self.quickprocess = kwargs.get('quickprocess',False)
		self.checkbad = kwargs.get('checkbad',False)
		self.meanVals = []
		self.rmsVals = []
		self.badVals = []
	def _preprocess(self,fits,f):
		self.imgMeans = []
		self.imgStds = []
		self.imgBad = []
	def process_hdu(self,extName,data,hdr):
		if self.checkbad:
			xy = np.indices(data.shape)
			inf = np.isinf(data)
			nan = np.isnan(data)
			data = np.ma.array(data,mask=np.logical_or(inf,nan))
			#zero = np.ma.less_equal(data,0)
			self.imgBad.append(xy[:,data.mask])
		if self.quickprocess:
			pix = overscan_subtract(data,hdr,method='mean_value',
			                        reject='sigma_clip',clip_iters=1,
			                        apply_filter=None)
		else:
			pix = data
		v,s = bokutil.array_stats(pix[self.statSec],method='mean',
		                          rms=True,clip=True,**self.clipArgs)
		self.imgMeans.append(v)
		self.imgStds.append(s)
		return data,hdr
	def _postprocess(self,fits,f):
		self.meanVals.append(self.imgMeans)
		self.rmsVals.append(self.imgStds)
		self.badVals.append(self.imgBad)
	def _finish(self):
		self.meanVals = np.array(self.meanVals)
		self.rmsVals = np.array(self.rmsVals)
	def reset(self):
		self.meanVals = []
		self.rmsVals = []
		self.badVals = []

def interpolate_masked_pixels(data,along='twod',method='linear'):
	if along=='rows':
		xx = np.arange(data.shape[1])
		for i in range(data.shape[0]):
			row = data.data[i]
			rowMask = data.mask[i]
			if rowMask.sum() > len(xx)-20:
				# less than 20 pixels on this row are not masked, don't
				# attempt interpolation
				row[rowMask] = np.nan
				continue
			interpFun = interp1d(xx[~rowMask],row[~rowMask],kind=method,
			                     bounds_error=False,fill_value=None)
			ii = np.where(rowMask)[0]
			rowInterp = interpFun(xx[ii])
			row[ii] = rowInterp
			## XXX fix extrapolation; for now don't overwrite those pix
			#good = np.where(~np.isnan(rowInterp))[0]
			#row[ii[good]] = rowInterp[good]
	elif along=='twod':
		y,x = np.indices(data.shape)
		yy = np.arange(data.shape[0])
		xx = np.arange(data.shape[1])
		interpFun = interp2d(x[~data.mask],y[~data.mask],
		                     data.data[~data.mask],kind='linear',
		                     bounds_error=False,fill_value=None)
		data.data[data.mask] = interpFun(x[data.mask],y[data.mask])
		print interpData.shape
		#interpData = griddata((y[~data.mask],x[~data.mask]),
		#                      data.data[~data.mask],
		#                      (yy,xx),method=method)
	else:
		raise ValueError
	return data

class BackgroundFit(object):
	def __init__(self,fits,nbin=64,coordsys='sky'):
		self.fits = fits
		self.binnedIm = fits.make_fov_image(nbin,coordsys,binfunc=np.ma.mean,
		                                    binclip=True,single=True,
		                                    mingood=nbin**2//3)
		self.coordSys = coordsys
	def __call__(self,extn):
		raise NotImplementedError
	def write(self,outFile,opfun=None,clobber=False):
		outFits = fitsio.FITS(outFile,'rw',clobber=clobber)
		outFits.write(None,header=self.fits.get_header(0))
		for extn,data,hdr in self.fits:
			im = self.get(extn)
			if opfun:
				im = opfun(im)
			outFits.write(im.astype(np.float32),extname=extn,header=hdr)
		outFits.close()

class SplineBackgroundFit(BackgroundFit):
	def __init__(self,fits,nKnots=2,order=1,**kwargs):
		super(SplineBackgroundFit,self).__init__(fits,**kwargs)
		self.nKnots = nKnots
		self.splineOrder = order
		im = self.binnedIm['im']
		x = self.binnedIm['x']
		y = self.binnedIm['y']
		tx = np.linspace(x.min(),x.max(),nKnots)
		ty = np.linspace(y.min(),y.max(),nKnots)
		self.spFit = LSQBivariateSpline(x[~im.mask],y[~im.mask],
		                                im.data[~im.mask],tx,ty,
		                                kx=self.splineOrder,
		                                ky=self.splineOrder)
	def _form_spline_im(self,ccdName,splinefun,xx,yy):
		# argh. spline requires monotonically increasing coordinates
		# as input, but my origin is at the center...
		if type(ccdName) is int:
			ccdName = 'CCD%d' % (ccdName+1)
		if ccdName == 'CCD1':
			ccdim = splinefun(xx[0,:],yy[:,0]).T
		elif ccdName == 'CCD2':
			ccdim = splinefun(xx[0,:],yy[::-1,0])[:,::-1].T
		elif ccdName == 'CCD3':
			ccdim = splinefun(xx[0,::-1],yy[:,0])[::-1,:].T
		elif ccdName == 'CCD4':
			ccdim = splinefun(xx[0,::-1],yy[::-1,0])[::-1,::-1].T
		return ccdim
	def __call__(self,x,y):
		return self.spFit(x,y)
	def get(self,extn):
		xi,yi = self.fits.get_xy(extn,self.coordSys)
		return self._form_spline_im(extn,self.spFit,xi,yi)

class PolynomialBackgroundFit(BackgroundFit):
	def __init__(self,fits,order=1,**kwargs):
		super(PolynomialBackgroundFit,self).__init__(fits,**kwargs)
		self.polyOrder = order
		im = self.binnedIm['im']
		x = self.binnedIm['x']
		y = self.binnedIm['y']
		self.polyModel = models.Polynomial2D(degree=order)
		self.polyFitFun = fitting.LinearLSQFitter()
		self.polyFit = self.polyFitFun(self.polyModel,
		                               x[~im.mask],y[~im.mask],
		                               im.data[~im.mask])
	def __call__(self,x,y):
		return self.polyFit(x,y)
	def get(self,extn):
		xi,yi = self.fits.get_xy(extn,self.coordSys)
		return self.polyFit(xi,yi)

class BokBiasStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','amp_central_quadrant')
		# this helps to deal with the roll-offs frequently seen in <=2015
		kwargs.setdefault('clip_sig',2.0)
		kwargs.setdefault('clip_cenfunc',np.ma.median)
		super(BokBiasStack,self).__init__(**kwargs)
		self.headerKey = 'BIAS'
		self.findRollOffs = kwargs.get('find_rolloffs',False)
		self.minNexp = 5
	def _reject_pixels(self,imCube):
		imCube = super(BokBiasStack,self)._reject_pixels(imCube)
		if self.findRollOffs:
			bslice = imCube[5:10,500:1500].reshape(-1,imCube.shape[-1])
			v = np.ma.median(bslice,axis=0).filled()
			bad = np.where(v < -50)[0]
			imCube[:100,:,bad] = np.ma.masked
		return imCube
	def _postprocess(self,extName,stack,hdr):
		colmedian = np.ma.median(stack,axis=0).filled(0)
		# don't know of a better way to fill along columns
		colmedian = np.repeat(colmedian[np.newaxis,:],stack.shape[0],axis=0)
		ismasked = stack.mask.copy() # copy suppresses a warning
		stack[ismasked] = colmedian[ismasked]
		return stack,hdr

class BokDomeFlatStack(bokutil.ClippedMeanStack):
	overExposedFlatCounts = 45000
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','amp_central_quadrant')
		kwargs.setdefault('scale','normalize_mode')
		super(BokDomeFlatStack,self).__init__(**kwargs)
		self.headerKey = 'FLAT'
		self.normPix = kwargs.get('norm_region',
		                    bokutil.stats_region('amp_corner_ccdcenter_256'))
	def _postprocess(self,extName,stack,hdr):
		flatNorm = bokutil.array_stats(stack[self.normPix])
		if flatNorm > self.overExposedFlatCounts:
			raise OverExposedFlatException("counts=%d"%flatNorm)
		stack /= flatNorm
		try:
			stack = stack.filled(1.0)
		except:
			pass
		hdr['FLATNORM'] = float(flatNorm)
		if self.scales is not None:
			for _i,_scl in enumerate(self.scales,start=1):
				hdr['FLTSCL%02d'%_i] = float(_scl)
		return stack,hdr

class NormalizeFlat(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','NORMFLT')
		super(NormalizeFlat,self).__init__(**kwargs)
		self.nbin = kwargs.get('nbin',16)
		self.nKnots = kwargs.get('nknots',2)
		self.splineOrder = kwargs.get('spline_order',3)
		self.flatFitMap = kwargs.get('_normed_flat_fit_map')
		self.binnedFlatMap = kwargs.get('_binned_flat_map')
		self.normedFlatFit = None
		self.binnedFlat = None
	def _preprocess(self,fits,f):
		if self.flatFitMap is not None:
			self.normedFlatFit = fitsio.FITS(self.flatFitMap(f),'rw')
			self.normedFlatFit.write(None,header=fits.get_header(0))
		if self.binnedFlatMap is not None:
			self.binnedFlat = fitsio.FITS(self.binnedFlatMap(f),'rw')
			self.binnedFlat.write(None,header=fits.get_header(0))
	def _postprocess(self,fits,f):
		if self.normedFlatFit is not None:
			self.normedFlatFit.close()
			self.normedFlatFit = None
		if self.binnedFlat is not None:
			self.binnedFlat.close()
			self.binnedFlat = None
	def process_hdu(self,extName,data,hdr):
		ny,nx = data.shape
		im = np.ma.masked_array(data,mask=((data<0.5) | (data>1.5)))
		# mask the edges
		margin = 15
		im.mask[:margin] = True
		im.mask[-margin:] = True
		im.mask[:,:margin] = True
		im.mask[:,-margin:] = True
		if extName == 'IM9':
			# this one has a long bad strip
			im.mask[:,:50] = True
		binnedIm = bokutil.rebin(im,self.nbin)
		binnedIm = bokutil.array_clip(binnedIm,axis=-1,
		                    clip_iters=3,clip_sig=2.2).mean(axis=-1)#.filled(1)
		x = np.arange(self.nbin/2,nx,self.nbin,dtype=np.float32)
		y = np.arange(self.nbin/2,ny,self.nbin,dtype=np.float32)
		xx,yy = np.meshgrid(x,y)
		tx = np.linspace(0,nx,self.nKnots)
		ty = np.linspace(0,ny,self.nKnots)
		spFit = LSQBivariateSpline(xx[~binnedIm.mask],yy[~binnedIm.mask],
		                           binnedIm.data[~binnedIm.mask],tx,ty,
		                           kx=self.splineOrder,ky=self.splineOrder)
		gradientIm = spFit(np.arange(nx),np.arange(ny)).T
		normedIm = im.data / gradientIm
		if self.normedFlatFit is not None:
			self.normedFlatFit.write(gradientIm.astype(np.float32),
			                         extname=extName,header=hdr)
		if self.binnedFlat is not None:
			self.binnedFlat.write(binnedIm.astype(np.float32),
			                      extname=extName,header=hdr)
		return normedIm.astype(np.float32),hdr

class BokGenerateDataQualityMasks(bokutil.BokProcess):
	# setting a relatively low value here because in 2016 a couple of the
	# amps started overflowing around this ADU level, but this could perhaps
	# be more refined
	satVal = 55000
	_procMsg = 'data quality mask for %s'
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','DQMASK')
		# this hack makes it possible to iterate over the 16-amp extensions
		# but only write out the 4-ccd mask at the end
		kwargs['read_only'] = True
		super(BokGenerateDataQualityMasks,self).__init__(**kwargs)
	def process_file(self,f):
		# need to overload process_file in order to check for output before
		# opening input image
		outf = self.outputNameMap(f)
		if os.path.exists(outf):
			if self.clobber:
				os.unlink(outf)
			else:
				self._proclog('data quality mask for %s already exists' % f)
				return None
		super(BokGenerateDataQualityMasks,self).process_file(f)
	def _preprocess(self,fits,f):
		super(BokGenerateDataQualityMasks,self)._preprocess(fits,f)
		self.hduData = []
	@staticmethod
	def _grow_saturated_blobs(ccdIm,saturated,minNsat=1000):
		yi,xi = np.indices(ccdIm.shape)
		# identify contiguous blogs associated with saturated pixels 
		# (i.e., stars) and count the number of saturated pixels in 
		# each blob, then sort by the largest number (brightest stars)
		satBlobs = measure.label(saturated,background=0)
		# create image with nominal mask
		for blob in range(1,satBlobs.max()+1):
			nSat = np.sum(satBlobs==blob)
			if nSat < minNsat:
				continue
			# a quick & hokey centering algorithm -- middle of the 
			# saturated blob!
			yextent = np.sum(satBlobs==blob,axis=0)
			# this gets messed up if the bleed trails extend to the end of
			# the image, leading to an extended pool near the edge
			xextent = np.sum(np.any(satBlobs==blob,axis=0))
			if xextent > 4000:
				yextent[:100] = 0
				yextent[-100:] = 0
			jj = np.where(yextent==yextent.max())[0]
			xc = xi[0,jj].mean()
			jj = np.where((satBlobs==blob)[:,int(xc)])[0]
			yc = yi[jj,int(xc)].mean()
			R = np.sqrt((xi-xc)**2 + (yi-yc)**2)
			# empirically found to be a reasonable minimum radius
			rad = pow(10,0.5*(np.log10(nSat)-np.log10(50)) + 1.6)
			ccdIm[R<rad] = np.ma.masked
		return ccdIm
	def process_hdu(self,extName,data,hdr):
		self.hduData.append(data)
		return data,hdr
	def _postprocess(self,fits,f):
		maskOut = fitsio.FITS(self.outputNameMap(f),'rw')
		maskOut.write(None,header=fits.get_header(0))
		for ccdNum,extGroup in enumerate(amp_iterator(),start=1):
			ims = [ self.hduData[ampNum-1] for ampNum in extGroup ]
			hdr = fits.get_header(bokCenterAmps[ccdNum-1])
			ccdIm,hdr = _orient_mosaic(hdr,ims,ccdNum,'center')
			saturated = (ccdIm > self.satVal) & ~ccdIm.mask
			maskIm = np.zeros(ccdIm.shape,dtype=np.int8)
			maskIm[ccdIm.mask] |= 1
			maskIm[saturated] |= 2
			# XXX this goes to gaincal? or remove it?
			# bright stars have swamped the image, mask it all
			if False: #saturated.sum() > 50000:
				ccdIm[:,:] = np.ma.masked
			else:
				ccdIm = self._grow_saturated_blobs(ccdIm,saturated)
			# negative value means "warning"
			maskIm[(maskIm==0)&ccdIm.mask] = -1
			maskOut.write(maskIm,extname='CCD%d'%ccdNum,header=hdr)
		maskOut.close()

class BokCCDProcess(bokutil.BokProcess):
	_procMsg = 'ccdproc %s'
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','CCDPROC')
		super(BokCCDProcess,self).__init__(**kwargs)
		self.fixPix = kwargs.get('fixpix',False)
		self.fixPixAlong = kwargs.get('fixpix_along','rows')
		self.fixPixMethod = kwargs.get('fixpix_method','linear')
		self.gainMultiply = kwargs.get('gain_multiply',True)
		self.inputGain = kwargs.get('input_gain',{ 'IM%d'%ampNum:g 
		                          for ampNum,g in zip(ampOrder,nominal_gain)})
		self.divideExpTime = kwargs.get('divide_exptime',False)
		# hacky way to say arithmetic is appropriate for inverse-var maps
		self.asWeight = kwargs.get('asweight',False)
		# store the maps to calibration images
		self.imTypes = ['bias','flat','ramp','illum','fringe','skyflat']
		self.calib = {}
		for imType in self.imTypes:
			cal = kwargs.get(imType)
			if cal is None:
				cal = bokdm.NullCalibrator()
			self.calib[imType] = cal
	def _preprocess(self,fits,f):
		super(BokCCDProcess,self)._preprocess(fits,f)
		hdrCards = {}
		for imType in self.imTypes:
			self.calib[imType].setTarget(f)
			calFn = self.calib[imType].getFileName()
			if calFn is not None:
				# e.g., hdr['BIASFILE'] = <filename>
				hdrKey = str(imType.upper()+'FILE')[:8]
				curPath = calFn
				if len(curPath) > 50:
					# trim the file path if it is long
					fn = ''
					while True:
						curPath,curEl = os.path.split(curPath)
						if len(fn) == 0:
							fn = curEl
						else:
							_fn = os.path.join(curEl,fn)
							if len(_fn) < 50:
								fn = _fn
							else:
								break
					fn = os.path.join('...',fn)
				else:
					fn = curPath
				hdrCards[hdrKey] = fn
		if self.divideExpTime:
			#hdrCards['BUNIT'] = 'counts/s'
			hdrCards['UNITS'] = 'counts/s'
			self.curExpTime = fits.outFits[0].read_header()['EXPTIME']
		fits.outFits[0].write_keys(hdrCards)
	def process_hdu(self,extName,data,hdr):
		# XXX hacky way to preserve arithmetic on masked pixels,
		#     but need to keep mask for fixpix
		if type(data) is np.ma.core.MaskedArray:
			data = data.data
		fscl = None
		bias = self.calib['bias'].getImage(extName)
		if bias is not None:
			data -= bias
		ramp = self.calib['ramp'].getImage(extName)
		if ramp is not None:
			data -= ramp
		for flatType in ['flat','illum','skyflat']:
			if flatType == 'skyflat':
				# ugly, but this is just where it goes
				fringe = self.calib['fringe'].getImage(extName)
				if fringe is not None:
					fscl = self.calib['fringe'].getFringeScale(extName,data)
					data -= fringe * fscl
			flat = self.calib[flatType].getImage(extName)
			if flat is not None:
				if self.asWeight:
					data *= flat**2  # inverse variance
				else:
					data /= flat     # image counts
		if self.fixPix:
			data = interpolate_masked_pixels(data,along=self.fixPixAlong,
			                                 method=self.fixPixMethod)
			# what should the fill value be?
			data = data.filled(0)
		if self.fixPix:
			hdr['FIXPIX'] = self.fixPixAlong
		if self.gainMultiply:
			data *= self.inputGain[extName]
			chNum = int(extName.replace('IM',''))
			hdr['GAIN'] = 1.0 # now in e-
			hdr['SATUR'] = saturation_dn * self.inputGain[extName]
			hdr['GAIN%02dA'%chNum] = self.inputGain[extName]
		elif 'SATUR' not in hdr:
			hdr['SATUR'] = saturation_dn
		if fscl is not None:
			hdr['FRNGSCL'] = float(fscl)
		if self.divideExpTime:
			data /= self.curExpTime
			hdr['SATUR'] /= self.curExpTime
			hdr['GAIN'] *= self.curExpTime
		return data,hdr

class BokWeightMap(bokutil.BokProcess):
	_procMsg = 'weight map %s'
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','WHTMAP')
		super(BokWeightMap,self).__init__(**kwargs)
		self._mask_map = kwargs.get('_mask_map')
		if isinstance(self._mask_map,fitsio.FITS):
			self._mask_map = bokutil.FakeFITS(self._mask_map)
		self.inputGain = kwargs.get('input_gain',{ 'IM%d'%ampNum:g 
		                          for ampNum,g in zip(ampOrder,nominal_gain)})
		self.maskFile = None
		self.flat = kwargs.get('flat')
		if self.flat is None:
			self.flat = bokdm.NullCalibrator()
	def _preprocess(self,fits,f):
		super(BokWeightMap,self)._preprocess(fits,f)
		try:
			maskFile = self._mask_map(f)
		except:
			maskFile = self._mask_map
		if maskFile != self.maskFile:
			if isinstance(maskFile,basestring):
				self.maskFile = self._mask_map(f)
				self.maskFits = fitsio.FITS(self.maskFile)
			else:
				self.maskFits = maskFile
		self.flat.setTarget(f)
		calFn = self.flat.getFileName()
		fits.outFits[0].write_keys({'FLATFILE':calFn})
	def process_hdu(self,extName,data,hdr):
		data,oscan_cols,oscan_rows = extract_overscan(data,hdr)
		data,mask = bokutil.mask_saturation(extName,data)
		mask |= ( (self.maskFits[extName][:,:] > 0) |
		          (data==0) )
		data -= np.median(oscan_cols)
#		if oscan_rows is not None:
#			data -= np.median(oscan_rows)
		flatField = self.flat.getImage(extName)
		if flatField is None:
			flatField = 1.0
		gain = self.inputGain[extName]
		ivar = np.clip(data,1e-10,65535)**-1
		ivar *= (gain/flatField)**-2
		ivar[mask] = 0
		#
		chNum = int(extName.replace('IM',''))
		hdr['GAIN%02dA'%chNum] = gain
		return ivar,hdr

class BokSkySubtract(bokutil.BokProcess):
	_procMsg = 'sky subtracting %s'
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','SKYSUB')
		super(BokSkySubtract,self).__init__(**kwargs)
		self.skyFitMap = kwargs.get('skyfit_map')
		self.skyFits = None
		self.method = kwargs.get('method','spline')
		self.order = kwargs.get('order',1)
		self.nBin = kwargs.get('nbin',64)
		self.nKnots = kwargs.get('nKnots',2)
		if self.method == 'spline':
			self.fitKwargs = {'nbin':self.nBin,'nKnots':self.nKnots,
			                  'order':self.order,'coordsys':'sky'}
			self.fitGen = SplineBackgroundFit
		elif self.method == 'polynomial':
			self.fitKwargs = {'nbin':self.nBin,
			                  'order':self.order,'coordsys':'sky'}
			self.fitGen = PolynomialBackgroundFit
		else:
			raise ValueError('fit method %s unknown' % method)
	def _fit_sky_model(self,fits):
		self.skyFit = self.fitGen(fits,**self.fitKwargs)
		# subtract the sky level at the origin so as to only remove 
		# the gradient
		self.sky0 = self.skyFit(0,0)
	def _preprocess(self,fits,f):
		super(BokSkySubtract,self)._preprocess(fits,f)
		self._fit_sky_model(fits)
		if self.skyFitMap is not None:
			self.skyFits = fitsio.FITS(self.skyFitMap(f),'rw')
			self.skyFits.write(None,header=fits.get_header(0))
		hdrCards = {'SKYMETH':self.method,'SKYORDR':self.order,
		            'SKYNBIN':self.nBin}
		if self.method=='spline':
			hdrCards['SKYNKNOT'] = self.nKnots
		fits.outFits[0].write_keys(hdrCards)
	def _postprocess(self,fits,f):
		if self.skyFits is not None:
			self.skyFits.close()
			self.skyFits = None
	def process_hdu(self,extName,data,hdr):
		if type(data) is np.ma.core.MaskedArray:
			data = data.data
		skyFit = (self.skyFit.get(extName) - self.sky0).astype(np.float32)
		data -= skyFit
		hdr['SKYVAL'] = float(self.sky0)
		# XXX should write fit parameters to header
		if self.skyFits is not None:
			self.skyFits.write(skyFit,extname=extName,header=hdr)
		return data,hdr

###############################################################################
#                                                                             #
#                               GAIN BALANCE                                  #
#            balances the amplifiers/CCDs with a gain correction              #
#                                                                             #
###############################################################################

class BokCalcGainBalanceFactors(bokutil.BokProcess):
	_procMsg = 'calculating gain balance factors for %s'
	maskFracThresh = 0.67
	ampMap = [ ((2,2,None),(2,0,'y'),(2,3,'x'),(0,1,'x')),
	           ((1,1,None),(1,0,'x'),(1,3,'y'),(3,2,'x')),
	           ((1,1,None),(1,0,'x'),(1,3,'y'),(3,2,'x')),
	           ((0,0,None),(0,1,'x'),(0,2,'y'),(2,3,'x')) ]
	# avoid correcting CCD3 from CCD1, because that is the glowing edge
	# the result is a CCW rotation starting from CCD1
	ccdMap = [(1,None,None,None),(2,2,7,'x'),(4,4,13,'y'),(3,14,11,'x')]
	nSplineKnots = 1
	splineOrder = 2
	nSplineRejIter = 1
	splineRejThresh = 2.5
	nFillEdge = 3
	def __init__(self,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokCalcGainBalanceFactors,self).__init__(**kwargs)
		self.bsMaskNameMap = kwargs.get('ccd_mask_map')
		self.bsMaskType = kwargs.get('mask_type','nonzero')
		self.statsMethod = kwargs.get('stats_method','mode')
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.clipArgs.setdefault('clip_iters',2)
		self.clipArgs.setdefault('clip_sig',2.5)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.mean)#np.ma.median)
		self.ratioClipArgs = {'clip_iters':2,'clip_sig':2.0}
		self.amplen = kwargs.get('amp_strip_length',1024)
		self.ampwid = kwargs.get('amp_strip_width',32)
		self.ampstride = kwargs.get('amp_strip_stride',8)
		self.ccdlen = kwargs.get('ccd_strip_length',1500)
		self.ccdwid = kwargs.get('ccd_strip_width',(10,50))
		self.ccdstride = kwargs.get('ccd_strip_stride',10)
		j1,j2 = self.ccdwid
		self.ampEdgeSlices = {
		  'x':np.s_[-self.amplen::self.ampstride,-self.ampwid:],
		  'y':np.s_[-self.ampwid:,-self.amplen::self.ampstride]
		}
		self.ccdEdgeSlices = {
		  'x':np.s_[-self.ccdlen::self.ccdstride,j1:j2],
		  'y':np.s_[j1:j2,-self.ccdlen::self.ccdstride]
		}
		# hugely downsample
		self.skyRegion = bokutil.stats_region('amp_central_quadrant',10)
		self.gainTrendMethod = kwargs.get('gain_trend_meth','spline')
		assert self.gainTrendMethod in ['median','spline']
		self.reset()
	def reset(self):
		self.files = []
		self.ampRelGains = []
		self.ccdRelGains = []
		self.allSkyVals = []
		self.allSkyRms = []
	def _preprocess(self,fits,f):
		super(BokCalcGainBalanceFactors,self)._preprocess(fits,f)
		if self.nProc > 1:
			# this prevents the return values from _getOutput from piling up
			# with duplicates when a subprocess is reused
			self.reset()
		self.files.append(f)
		self.hduData = []
		self.rawSky = []
		self.rawSkyRms = []
	def process_hdu(self,extName,data,hdr):
		self.hduData.append(data)
		sky = bokutil.array_clip(data[self.skyRegion],clip_iters=2)
		self.rawSky.append(sky.mean())
		self.rawSkyRms.append(sky.std())
		return data,hdr
	def process_files(self,files,filters):
		self.filters = filters
		super(BokCalcGainBalanceFactors,self).process_files(files)
	def _slice_ok(self,imslice,nslice):
		return np.sum(imslice.mask)/float(nslice) < self.maskFracThresh
	def _get_img_slice(self,imslice):
		nslice = imslice.size
		if not self._slice_ok(imslice,nslice): 
			return None
		# XXX use the global sky+rms here instead of clipping within slice
		imslice = bokutil.array_clip(imslice,**self.clipArgs)
		imslice.mask |= binary_dilation(imslice.mask,iterations=1)
		if not self._slice_ok(imslice,nslice): 
			return None
		return imslice
	def _get_img_slices(self,refExt,calExt,edgedir,which):
		if which=='amp':
			s = self.ampEdgeSlices[edgedir]
		elif which=='ccd':
			s = self.ccdEdgeSlices[edgedir]
		refslice = self._get_img_slice(self.hduData[refExt][s])
		calslice = self._get_img_slice(self.hduData[calExt][s])
		if refslice is None or calslice is None:
			print 'could not correct ',which,refExt,calExt
		return refslice,calslice
	def _get_flux_ratio(self,refExt,calExt,edgedir,which):
		refslice,calslice = self._get_img_slices(refExt,calExt,edgedir,which)
		if refslice is None or calslice is None:
			return 0.0
		axis = {'x':1,'y':0}[edgedir]
		refv = refslice.mean(axis=axis)
		calv = calslice.mean(axis=axis)
		fratio = np.ma.divide(refv,calv)
		return bokutil.array_clip(fratio,**self.ratioClipArgs).mean()
	def _apply_bright_star_mask(self,f):
		bsMask = fitsio.FITS(self.bsMaskNameMap(f))
		for ccdNum,extGroup in enumerate(amp_iterator(),start=1):
			ccdMaskIm = bokutil.load_mask(bsMask['CCD%d'%ccdNum].read(),
			                              self.bsMaskType)
			ampMasks = bokutil.ccd_split(ccdMaskIm,ccdNum)
			for ampNum,ampMask in zip(extGroup,ampMasks):
				self.hduData[ampNum-1].mask |= ampMask
	def _postprocess(self,fits,f):
		self._apply_bright_star_mask(f)
		gains = np.zeros(16,dtype=np.float32)
		# balance the amps
		for ccdi,extGroup in enumerate(amp_iterator()):
			for ampi,ampj,edgedir in self.ampMap[ccdi]:
				refExt = 4*ccdi + ampi
				calExt = 4*ccdi + ampj
				if refExt==calExt:
					# the reference amp
					gains[calExt] = 1.0
				else:
					gains[calExt] = self._get_flux_ratio(refExt,calExt,
					                                     edgedir,'amp')
		# store the CCD count ratios
		ccdratios = np.ones(4,dtype=np.float32)
		for ccdNum,refExt,calExt,edgedir in self.ccdMap:
			if refExt==calExt:
				# the reference CCD
				ccdratios[ccdNum-1] = 1.0
			else:
				ccdratios[ccdNum-1] = self._get_flux_ratio(refExt,calExt,
				                                           edgedir,'ccd')
		#
		self.ampRelGains.append(gains)
		self.ccdRelGains.append(ccdratios)
		self.allSkyVals.append(np.array(self.rawSky))
		self.allSkyRms.append(np.array(self.rawSkyRms))
	def _getOutput(self):
		return (self.files,self.ampRelGains,self.ccdRelGains,self.allSkyVals)
	def _null_result(self,f):
		return ([f],[np.zeros(16,dtype=np.float32)],
		        [np.zeros(4,dtype=np.float32)],
		        [np.zeros(16,dtype=np.float32)])
	def _ingestOutput(self,procOut):
		self.files,self.ampRelGains,self.ccdRelGains,self.allSkyVals = \
		              zip(*procOut)
		# squeeze out the extra axis that occurs since output elements
		# are in lists of unit length when multiprocessing
		self.files = np.array(self.files).squeeze()
		self.filters = np.array(self.filters).squeeze()
		self.ampRelGains = np.array(self.ampRelGains).squeeze()
		self.ccdRelGains = np.array(self.ccdRelGains).squeeze()
		self.allSkyVals = np.array(self.allSkyVals).squeeze()
	def _median_gain_trend(self,gc):
		gc = np.ma.array(gc,mask=(gc==0),copy=True)
		gc_clip = sigma_clip(gc,iters=2,sigma=2.0,axis=0)
		gc[:] = gc_clip.mean(axis=0)
		return gc.filled(0),gc.mask
	def _spline_gain_trend(self,rawgc):
		nimg,ngain = rawgc.shape
		seqno = np.arange(nimg,dtype=np.float32)
		gc = np.ma.array(rawgc,mask=(rawgc==0),copy=True)
		msk = gc.mask.copy()
		knots = np.linspace(0,nimg,self.nSplineKnots+2)[1:-1]
		for j in range(ngain):
			if np.allclose(gc[:,j],1):
				# a reference chip
				continue
			ii = np.where(~gc[:,j].mask)[0]
			rejmask = np.zeros(len(seqno),dtype=np.bool)
			try:
				spfit = LSQUnivariateSpline(seqno[ii],gc[ii,j].filled(),
				                            knots,bbox=[0,nimg],
				                            k=self.splineOrder)
				for iternum in range(self.nSplineRejIter):
					ii = np.where(~(gc[:,j].mask | rejmask))[0]
					resid = gc[ii,j] - spfit(seqno[ii])
					resrms = np.ma.std(resid)
					rejmask[ii] |= np.abs(resid/resrms) > self.splineRejThresh
					ii = np.where(~(gc[:,j].mask | rejmask))[0]
					spfit = LSQUnivariateSpline(seqno[ii],gc[ii,j].filled(),
					                            knots,bbox=[0,nimg],
					                            k=self.splineOrder)
				gc[:,j] = spfit(seqno)
				msk[:,j] |= rejmask
				if ii[0] >= self.nFillEdge:
					# force it to linear fit with no interior knots
					i0 = min(5,len(ii)-5)
					linspfit = LSQUnivariateSpline(seqno[ii[:i0]],
					                               gc[ii[:i0],j].filled(),
					                               [],bbox=[0,seqno[ii[i0]]],
					                               k=1)
					gc[:ii[0],j] = linspfit(seqno[:ii[0]])
				if ii[-1] < len(seqno)-self.nFillEdge:
					i0 = max(0,len(ii)-5)
					linspfit = LSQUnivariateSpline(seqno[ii[i0:]],
					                               gc[ii[i0:],j].filled(),
					                               [],bbox=[seqno[ii[i0]],
					                                        seqno[ii[-1]]],
					                               k=1)
					gc[ii[-1]+1:,j] = linspfit(seqno[ii[-1]+1:])
			except ValueError:
				print 'WARNING: spline fit failed, reverting to mean'
				gc[:,j] = sigma_clip(gc[:,j],iters=2,sigma=2.0).mean()
		return gc.filled(0),msk
	def calc_mean_corrections(self):
		raw_ampg = self.ampRelGains = np.array(self.ampRelGains)
		raw_ccdg = self.ccdRelGains = np.array(self.ccdRelGains)
		filts = np.unique(self.filters)
		ampg = np.zeros_like(raw_ampg)
		ccdg = np.zeros_like(raw_ccdg)
		for filt in filts:
			ii = np.where(self.filters == filt)[0]
			if len(ii)==0: continue
			if self.gainTrendMethod == 'median':
				ampg[ii],msk = self._median_gain_trend(raw_ampg[ii])
				ccdg[ii],msk = self._median_gain_trend(raw_ccdg[ii])
			elif self.gainTrendMethod == 'spline':
				ampg[ii],msk = self._spline_gain_trend(raw_ampg[ii])
				ccdg[ii],msk = self._spline_gain_trend(raw_ccdg[ii])
		# propagate the gain corrections starting from the reference
		ampgscale = ampg.copy()
		for ccdi,extGroup in enumerate(amp_iterator()):
			for ampi,ampj,edgedir in self.ampMap[ccdi]:
				refExt = 4*ccdi + ampi
				calExt = 4*ccdi + ampj
				ampgscale[:,calExt] *= ampgscale[:,refExt]
		ccdgscale = ccdg.copy()
		for ccdNum,refExt,calExt,edgedir in self.ccdMap:
			if refExt!=calExt:
				refCcd = refExt // 4
				ccdgscale[:,ccdNum-1] *= ccdgscale[:,refCcd] * \
				                    (ampgscale[:,refExt]/ampgscale[:,calExt]) 
		self.ampGainTrend = ampg
		self.ccdGainTrend = ccdg
		self.gainCors = np.dstack([ampgscale,
		                           np.repeat(ccdgscale,4,axis=1)])
		self.correctedGains = np.dstack([raw_ampg*(ampgscale/ampg),
		                        np.repeat(raw_ccdg*(ccdgscale/ccdg),4,axis=1)])
		return self.gainCors
	def get_values(self):
		return ( np.array(self.ampRelGains),
		         np.array(self.ccdRelGains),
		         np.array(self.correctedGains),
		         np.array(self.ampGainTrend),
		         np.array(self.ccdGainTrend),
		         np.array(self.allSkyVals) )

###############################################################################
#                                                                             #
#                               COMBINE CCDs                                  #
#                balances the amplifiers with a gain correction               #
#                                                                             #
###############################################################################

def _orient_mosaic(hdr,ims,ccdNum,origin):
	outIm = bokutil.ccd_join(ims,ccdNum,origin=origin)
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
		hdr['CTYPE1'] = 'RA---TAN'
		hdr['CTYPE2'] = 'DEC--TAN'
		crval1 = hdr['CRVAL1'] 
		hdr['CRVAL1'] = hdr['CRVAL2']
		hdr['CRVAL2'] = crval1
		signx = [+1,-1][det_i]
		signy = [-1,+1][det_j]
		hdr['CD1_2'] = 0.0
		hdr['CD2_1'] = 0.0
		hdr['CD1_1'] = signx*plateScale
		hdr['CD2_2'] = signy*plateScale
		hdr['CRPIX1'] = -182.01
		hdr['CRPIX2'] = -59.04
		hdr['LTM1_1'] = -float(signx)
		hdr['LTM2_2'] = float(signy)
		hdr['LTV1'] = -182.
		hdr['LTV2'] = -59.
	elif origin == 'lower left':
		raise NotImplementedError # this is vestigial, need to update
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

# XXX should fit this into a BokProcess even if it requires some mungeing
#     of the way extensions are combined into ccds

def _combine_ccds(f,**kwargs):
	inputFileMap = kwargs.get('input_map',IdentityNameMap)
	outputFileMap = kwargs.get('output_map',IdentityNameMap)
	gainMap = kwargs.get('gain_map')
	origin = kwargs.get('origin','center')
	clobber = kwargs.get('clobber')
	# a hacked entry point for flat fields that normalizes the corners
	# of each channel at the CCD centers to be unity
	flatNorm = kwargs.get('apply_flat_norm',False)
	ignoreExisting = kwargs.get('ignore_existing',True)
	# another hacky entry point for preprocessing of per-amp images
	# before combining
	_preprocess_ims = kwargs.get('_preprocess_function')
	# fitsio doesn't accept file descriptors, but tempfile opens files...
	tmpFile = tempfile.NamedTemporaryFile()
	tmpFileName = tmpFile.name
	tmpFile.close()
	# do the extensions in numerical order, instead of HDU list order
	extns = np.array(['IM%d' % ampNum for ampNum in range(1,17)])
	#
	inputFile = inputFileMap(f)
	outputFile = outputFileMap(f)
	bokutil.mplog('combine_ccds: '+f,kwargs.get('processes',1))
	inFits = fitsio.FITS(inputFile)
	if 'CCDJOIN' in inFits[0].read_header():
		print '%s already combined, skipping' % f
		inFits.close()
		return
	if outputFile != inputFile:
		if os.path.exists(outputFile):
			if clobber:
				os.unlink(outputFile)
			# XXX should be checking header key here?
			elif ignoreExisting:
				print '%s already exists, skipping' % outputFile
				return
			else:
				raise bokutil.OutputExistsError(
				                    '%s already exists'%outputFile)
		outFits = fitsio.FITS(outputFile,'rw')
	else:
		# have to use a temporary file to change format
		if os.path.exists(tmpFileName):
			os.unlink(tmpFileName)
		outFits = fitsio.FITS(tmpFileName,'rw')
	hdr = inFits[0].read_header()
	hdr['DETSIZE'] = '[1:%d,1:%d]' % (8192,8064) # hardcoded
	hdr['NEXTEND'] = 4
	hdr['CCDJOIN'] = bokutil.get_timestamp()
	outFits.write(None,header=hdr)
	refSkyCounts = None
	for ccdNum,extGroup in enumerate(np.hsplit(extns,4),start=1):
		hdr = inFits[bokCenterAmps[ccdNum-1]].read_header()
		ccdIms = []
		satvals = []
		for j,ext in enumerate(extGroup):
			im = inFits[ext].read() 
			hext = inFits[ext].read_header()
			if _preprocess_ims is not None:
				im = _preprocess_ims(im,hext,ext)
			try:
				# copy in the nominal gain values from per-amp headers
				gainKey = 'GAIN%02dA' % int(ext[2:])
				g0 = hext[gainKey]
				hdr[gainKey] = g0
			except:
				pass
			if gainMap is not None:
				ampIdx = ampOrder[4*(ccdNum-1)+j] - 1
				gc1,gc2 = gainMap['corrections'][f][ampIdx]
				sky = gainMap['skyvals'][f][ampIdx]
				hdr['SKY%02dB'%int(ext[2:])] = sky
				hdr['GAIN%02dB'%int(ext[2:])] = gc1
				if j==0:
					hdr['CCDGAIN'] = gc2
					try:
						satvals.append(hdr['SATUR']*gc1*gc2)
					except ValueError:
						pass
				im *= gc1 * gc2
				# the final gain factor
				hdr['GAIN%02d'%int(ext[2:])] = g0 * gc1 * gc2
			if flatNorm:
				_s = bokutil.stats_region('amp_corner_ccdcenter_128')
				_a = bokutil.array_stats(im[_s],method='mode')
				im /= _a
			ccdIms.append(im)
		# orient the channel images into a mosaic of CCDs and
		# modify WCS & mosaic keywords
		outIm,hdr = _orient_mosaic(hdr,ccdIms,ccdNum,origin)
		if len(satvals)>0:
			hdr['SATUR'] = np.min(satvals)
		if True:
			# For some reason this results in a segfault when running
			# on NERSC unless the image is copied...
			outIm = outIm.copy()
		outFits.write(outIm,extname='CCD%d'%ccdNum,header=hdr)
	outFits.close()
	if outputFile == inputFile:
		shutil.move(tmpFileName,inputFile)

def _combine_ccds_exc(f,**kwargs):
	try:
		_combine_ccds(f,**kwargs)
	except Exception,e:
		sys.stderr.write('FAILED: combine %s [%s]\n'%(f,e))

def combine_ccds(fileList,**kwargs):
	procmap = kwargs.pop('procmap',map)
	if kwargs.get('debug',False):
		combfunc = partial(_combine_ccds,**kwargs)
	else:
		combfunc = partial(_combine_ccds_exc,**kwargs)
	procmap(combfunc,fileList)


###############################################################################
#                                                                             #
#                              PROCESS ROUND 2                                #
#   obj detection & masking, divide by supersky flat, ... (CR rejection?)     #
#                                                                             #
###############################################################################

def find_bright_stars(im,saturation,minNSat=100):
	y,x = np.indices(im.shape)
	saturated = im >= saturation
	satObjs,nObjs = meas.label(saturated)#,np.ones((3,3)))
	nsat = meas.labeled_comprehension(saturated,satObjs,
	                                  np.arange(1,nObjs+1),np.sum,int,-1)
	ii = np.where(nsat > minNSat)[0]
	# no way to only measure by index?
	blobs = meas.find_objects(satObjs)
	cntrx = [ np.average(x[blobs[i]][~saturated[blobs[i]]],
	                     weights=im[blobs[i]][~saturated[blobs[i]]])
	            for i in ii ]
	cntry = [ np.average(y[blobs[i]][~saturated[blobs[i]]],
	                     weights=im[blobs[i]][~saturated[blobs[i]]])
	            for i in ii ]
	return np.array([cntrx,cntry]).astype(int)

def mask_bright_stars(im,saturation,minNSat=50):
	y,x = np.indices(im.shape)
	mask = np.zeros(im.shape,dtype=bool)
	saturated = im >= saturation
	satObjs,nObjs = meas.label(saturated)#,np.ones((3,3)))
	nsat = meas.labeled_comprehension(saturated,satObjs,
	                                  np.arange(1,nObjs+1),np.sum,int,-1)
	ii = np.where(nsat > minNSat)[0]
	# no way to only measure by index?
	blobs = meas.find_objects(satObjs)
	for i in ii:
		cntrx = np.average(x[blobs[i]][~saturated[blobs[i]]],
		                   weights=im[blobs[i]][~saturated[blobs[i]]])
		cntry = np.average(y[blobs[i]][~saturated[blobs[i]]],
		                   weights=im[blobs[i]][~saturated[blobs[i]]])
		cntrx = int(cntrx)
		cntry = int(cntry)
		xextent = np.max(np.abs(x[blobs[i]]-cntrx))
		mask[(cntry-xextent):(cntry+xextent),
		     (cntrx-xextent):(cntrx+xextent)] = True
	return mask

class BokGenerateSkyFlatMasks(bokutil.BokProcess):
	_procMsg = 'generating sky mask for %s'
	def __init__(self,**kwargs):
		self.nBin = kwargs.get('binSize',4)
		kwargs.setdefault('header_key','SKYMSK')
		super(BokGenerateSkyFlatMasks,self).__init__(**kwargs)
		self.hiThresh = kwargs.get('high_thresh',5.0)
		self.loThresh = kwargs.get('low_thresh',5.0)
		self.growThresh = kwargs.get('grow_thresh',1.0)
		self.binGrowSize = kwargs.get('mask_grow_size',3)
		self.nSample = kwargs.get('num_sample',16)
		self.nKnots = kwargs.get('num_spline_knots',3)
		self.splineOrder = kwargs.get('spline_order',3)
		self.statsPix = bokutil.stats_region(kwargs.get('stats_region'),
		                                     self.nSample)
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.clipArgs.setdefault('clip_iters',3)
		self.clipArgs.setdefault('clip_sig',2.2)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.mean)
		self.growKern = None #np.ones((self.binGrowSize,self.binGrowSize),dtype=bool)
		self.nPad = 10
		self.noConvert = True
	def process_hdu(self,extName,data,hdr):
		if (data>hdr['SATUR']).sum() > 50000:
			# if too many pixels are saturated mask the whole damn thing
			hdr['BADSKY'] = 1
		sky,rms = bokutil.array_stats(data[self.statsPix],
		                              method='mode',rms=True,
		                              **self.clipArgs)
		binnedIm = bokutil.rebin(data,self.nBin)
		# propagate the mask if too many sub-pixels are masked
		#mask = binnedIm.sum(axis=-1) > self.nBin**2/2
		binnedIm = binnedIm.mean(axis=-1)
		mask = binnedIm.mask.copy()
		# divide by RMS to make a SNR image
		snr = (binnedIm.data-sky) / (rms/self.nBin)
		# fit and subtract a smooth sky model, after a first round of
		# masking sources 
		mask |= (snr < -15.0) | (snr > 15.0) 
		# construct a spline model for the sky background
		y,x = np.indices(binnedIm.shape)
		tx = np.linspace(0,binnedIm.shape[1],self.nKnots)
		ty = np.linspace(0,binnedIm.shape[0],self.nKnots)
		spfit = LSQBivariateSpline(x[~mask],y[~mask],binnedIm[~mask],
		                       tx,ty,kx=self.splineOrder,ky=self.splineOrder)
		# remake the SNR image after subtracting fitted sky
		snr = (binnedIm.data - spfit(x[0],y[:,0]).T) / (rms/self.nBin)
		# mask everything above the threshold
		mask |= (snr < -self.loThresh) | (snr > self.hiThresh)
		# grow the mask from positive deviations
		binary_dilation(mask,mask=(snr>self.growThresh),iterations=0,
		                structure=self.growKern,output=mask)
		# fill in holes in the mask -- binary closing will lose masked
		# values at edge, so have to pad the mask first
		# should also check scipy.ndimage.morphology.binary_fill_holes()...
		maskpad = np.pad(mask,self.nPad,mode='constant',constant_values=0)
		binary_closing(maskpad,iterations=5,structure=self.growKern,
		               output=maskpad)
		mask[:] |= maskpad[self.nPad:-self.nPad,self.nPad:-self.nPad]
		# bright star mask
		bsmask = False # XXX
		# construct the output array
		maskIm = bokutil.magnify(mask,self.nBin) | bsmask
		return maskIm.astype(np.uint8),hdr

class BokFringePatternStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','ccd_central_quadrant')
		self.nSample = kwargs.get('num_sample',4)
		#kwargs.setdefault('scale','normalize_mode')
		kwargs.setdefault('maxmem',5)
		kwargs.setdefault('fill_value',1.0)
		super(BokFringePatternStack,self).__init__(**kwargs)
		self.statsPix = bokutil.stats_region(kwargs.get('stats_region'),
		                                     self.nSample)
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.statsMethod = kwargs.get('stats_method','mode')
		self.clipArgs.setdefault('clip_iters',3)
		self.clipArgs.setdefault('clip_sig',2.2)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.mean)
		self.smoothingLength = kwargs.get('smoothing_length',0.05)
		self.rawStackFile = kwargs.get('raw_stack_file')
		self.rawStackFits = None
		self.procmap = kwargs.get('procmap',map)
		self.headerKey = 'FRG'
	def _getnorm(self,f):
		fits = bokutil.BokMefImage(self.inputNameMap(f),
		                           mask_file=self.maskNameMap(f),
		                           read_only=True)
		meanVals = []
		for extn,data,hdr in fits:
			meanVal = bokutil.array_stats(data[self.statsPix],
			                              method=self.statsMethod,
			                              **self.clipArgs)
			meanVals.append(meanVal)
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] '%pid,
		print 'medians for image %s are %.1f %.1f %.1f %.1f' % \
		           tuple((self.inputNameMap(f),)+tuple(meanVals))
		return meanVals
	def _preprocess(self,fileList,outFits):
		# calculate the norms in subprocesses, need the workaround to
		# avoid sending a Pool object
		procmap = self.procmap
		self.procmap = None
		norms = procmap(self._getnorm,fileList)
		self.norms = np.array(norms).astype(np.float32)
		self.procmap = procmap
		if self.rawStackFile is not None:
			print 'writing raw stack to ',self.rawStackFile(outFits._filename)
			rawFn = self.rawStackFile(outFits._filename)
			if os.path.exists(rawFn):
				os.remove(rawFn)
			self.rawStackFits = fitsio.FITS(rawFn,'rw')
			# if we've gotten to here, we already know any existing file 
			# needs to be clobbered (XXX but this didn't work??? added above)
			self.rawStackFits.write(None,header=outFits[0].read_header(),
			                        clobber=True)
	def _rescale(self,imCube,scales=None):
		if scales is not None:
			_scales = scales[np.newaxis,:]
		else:
			_scales = self.norms[np.newaxis,:]
		_scales = _scales.mean(axis=-1) # XXX averaging across CCDs
		self.scales = _scales.squeeze()
		return imCube - _scales
	def _postprocess(self,extName,stack,hdr):
		cleanStack = stack.filled(1.0)
		interpMask = False
		if self.rawStackFile is not None:
			self.rawStackFits.write(cleanStack,extname=extName,header=hdr)
		cleanStack = spline_filter(cleanStack,self.smoothingLength)
		# renormalize to unity, using the combined interp and input mask
		_stack = np.ma.masked_array(cleanStack,mask=interpMask|stack.mask)
		normpix = bokutil.array_clip(_stack[self.statsPix])
		_stack -= normpix.mean()
		return _stack,hdr
	def _cleanup(self):
		super(BokFringePatternStack,self)._cleanup()
		if self.rawStackFits is not None:
			self.rawStackFits.close()
			self.rawStackFits = None

class BokNightSkyFlatStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','ccd_central_quadrant')
		kwargs.setdefault('scale','normalize_mode')
		kwargs.setdefault('maxmem',5)
		kwargs.setdefault('fill_value',1.0)
		super(BokNightSkyFlatStack,self).__init__(**kwargs)
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		# override some defaults
		self.statsPix = bokutil.stats_region(self.statsRegion,8)
		self.statsMethod = kwargs.get('stats_method','mode')
		self.clipArgs['clip_iters'] = 3
		self.clipArgs['clip_sig'] = 2.2
		self.clipArgs['clip_cenfunc'] = np.ma.mean
		self.smoothingLength = kwargs.get('smoothing_length',0.05)
		self.procmap = kwargs.get('procmap',map)
		self.normCCD = 'CCD1'
		self.headerKey = 'SKYFL'
	def _getnorm(self,f):
		fits = bokutil.BokMefImage(self.inputNameMap(f),
		                           mask_file=self.maskNameMap(f),
		                           read_only=True)
		# XXX have the get the full image and then subsample, because
		#     fitsio doesn't handle negative slice boundaries
		normpix = fits.get(self.normCCD)[self.statsPix]
		meanVal = bokutil.array_stats(normpix,method=self.statsMethod,
		                              **self.clipArgs)
		norm = 1/meanVal
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] '%pid,
		print 'norm for image %s is %f' % \
		           (self.inputNameMap(f),meanVal)
		return norm
	def _preprocess(self,fileList,outFits):
		# calculate the norms in subprocesses, need the workaround to
		# avoid sending a Pool object
		procmap = self.procmap
		self.procmap = None
		norms = procmap(self._getnorm,fileList)
		self.norms = np.array(norms).astype(np.float32)
		self.procmap = procmap
	def _rescale(self,imCube,scales=None):
		if scales is not None:
			_scales = scales[np.newaxis,:]
		else:
			_scales = self.norms[np.newaxis,:]
		self.scales = _scales.squeeze()
		return imCube * _scales
	def _postprocess(self,extName,stack,hdr):
		# renormalize to unity
		stack /= bokutil.array_clip(stack[self.statsPix]).mean()
		return stack,hdr



