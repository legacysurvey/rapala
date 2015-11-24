#!/usr/bin/env python

import os
import re
import shutil
import tempfile
import subprocess
import numpy as np
from scipy.interpolate import LSQBivariateSpline,RectBivariateSpline,griddata
from scipy.signal import spline_filter
from astropy.stats import sigma_clip
from astropy.modeling import models,fitting
import fitsio

from . import bokutil

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
    1.42733335,  1.38764536,  1.3538, 1.45403028
  ] )

saturation_dn = 65000

# XXX
configdir = os.environ['BOKPIPE']+'/config/'

from scipy.interpolate import interp1d,interp2d

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
		extensions = ['CCD%d'%i for i in range(1,5)] # XXX 4-CCD only...
	def __call__(self,extn):
		raise NotImplementedError

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
		super(BokBiasStack,self).__init__(**kwargs)
		self.headerKey = 'BIAS'

class BokDomeFlatStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','amp_corner_ccdcenter_small')
		kwargs.setdefault('scale','normalize_mode')
		super(BokDomeFlatStack,self).__init__(**kwargs)
		self.headerKey = 'FLAT'
	def _postprocess(self,extName,stack,hdr):
		flatNorm = bokutil.array_stats(stack[self.statsPix])
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
		self.nbin = kwargs.get('nbin',32)
		self.flatFitName = kwargs.get('normed_flat_fit_file')
		self.binnedFlatName = kwargs.get('binned_flat_file')
		self.normedFlatFit = None
		self.binnedFlat = None
	def _preprocess(self,fits,f):
		if self.flatFitName is not None:
			if os.path.exists(self.flatFitName):
				os.unlink(self.flatFitName)
			self.normedFlatFit = fitsio.FITS(self.flatFitName,'rw')
			self.normedFlatFit.write(None,header=fits.get_header(0))
		if self.binnedFlatName is not None:
			if os.path.exists(self.binnedFlatName):
				os.unlink(self.binnedFlatName)
			self.binnedFlat = fitsio.FITS(self.binnedFlatName,'rw')
			self.binnedFlat.write(None,header=fits.get_header(0))
	def _finish(self):
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
		binnedIm = sigma_clip(binnedIm,axis=-1,iters=3,sig=2.2,
		                      cenfunc=np.ma.mean).mean(axis=-1).filled(1)
		x = np.arange(self.nbin/2,nx,self.nbin)
		y = np.arange(self.nbin/2,ny,self.nbin)
		spfit = RectBivariateSpline(x,y,binnedIm.T,s=1)
		gradientIm = spfit(np.arange(nx),np.arange(ny)).T
		normedIm = im.data / gradientIm
		if self.normedFlatFit is not None:
			self.normedFlatFit.write(gradientIm.astype(np.float32),
			                         extname=extName,header=hdr)
		if self.binnedFlat is not None:
			self.binnedFlat.write(binnedIm.astype(np.float32),
			                      extname=extName,header=hdr)
		return normedIm.astype(np.float32),hdr

class BokCCDProcess(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','CCDPROC')
		super(BokCCDProcess,self).__init__(**kwargs)
		self.fixPix = kwargs.get('fixpix',False)
		self.fixPixAlong = kwargs.get('fixpix_along','rows')
		self.fixPixMethod = kwargs.get('fixpix_method','linear')
		self.gainMultiply = kwargs.get('gain_multiply',True)
		self.inputGain = kwargs.get('input_gain',{ 'IM%d'%ampNum:g 
		                          for ampNum,g in zip(ampOrder,nominal_gain)})
		self.imTypes = ['bias','flat','ramp','illum','darksky']
		#
		self.procIms = {}
		for imType in self.imTypes:
			fitsIn = kwargs.get(imType)
			self.procIms[imType] = {'file':None,'fits':None}
			self.procIms[imType]['master'] = True
			if type(fitsIn) is fitsio.fitslib.FITS:
				self.procIms[imType]['file'] = fitsIn._filename
				self.procIms[imType]['fits'] = fitsIn
			elif type(fitsIn) is str:
				self.procIms[imType]['file'] = fitsIn
				self.procIms[imType]['fits'] = fitsio.FITS(fitsIn)
			elif type(fitsIn) is dict:
				self.procIms[imType]['map'] = fitsIn
				self.procIms[imType]['master'] = False
	def _preprocess(self,fits,f):
		print 'ccdproc ',fits.fileName,fits.outFileName
		for imType in self.imTypes:
			if not self.procIms[imType]['master']:
				inFile = self.procIms[imType]['map'][f]
				if self.procIms[imType]['file'] != inFile:
					if self.procIms[imType]['fits'] is not None:
						self.procIms[imType]['fits'].close()
					self.procIms[imType]['file'] = inFile
					self.procIms[imType]['fits'] = fitsio.FITS(inFile)
			hdrCards = {}
			if self.procIms[imType]['file'] is not None:
				# e.g., hdr['BIASFILE'] = <filename>
				hdrKey = str(imType.upper()+'FILE')[:8]
				curPath = self.procIms[imType]['file'].rstrip('.fits')
				if len(curPath) > 65:
					# trim the file path if it is long
					fn = ''
					while len(fn)<60:
						curPath,curEl = os.path.split(curPath)
						if len(fn) == 0:
							fn = curEl
						else:
							fn = os.path.join(curEl,fn)
					fn = os.path.join('...',fn)
				else:
					fn = curPath
				hdrCards[hdrKey] = fn
			fits.outFits[0].write_keys(hdrCards)
	def process_hdu(self,extName,data,hdr):
		bias = self.procIms['bias']['fits'] 
		if bias is not None:
			data -= bias[extName][:,:]
		ramp = self.procIms['ramp']['fits'] 
		if ramp is not None:
			# this correction may not exist on all extensions
			try:
				data -= ramp[extName][:,:]
			except:
				pass
		for flatType in ['flat','illum','darksky']:
			flat = self.procIms[flatType]['fits']
			if flat is not None:
				data /= flat[extName][:,:]
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
		else:
			hdr['SATUR'] = saturation_dn
		return data,hdr

class BokSkySubtract(bokutil.BokProcess):
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
		print 'sky subtracting ',fits.fileName,fits.outFileName
		self._fit_sky_model(fits)
		if self.skyFitMap is not None:
			self.skyFits = fitsio.FITS(self.skyFitMap(f),'rw')
			self.skyFits.write(None,header=fits.get_header(0))
		hdrCards = {'SKYMETH':self.method,'SKYORDR':self.order,
		            'SKYNBIN':self.nBin}
		if self.method=='spline':
			hdrCards['SKYNKNOT'] = self.nKnots
		fits.outFits[0].write_keys(hdrCards)
	def process_hdu(self,extName,data,hdr):
		skyFit = (self.skyFit.get(extName) - self.sky0).astype(np.float32)
		data -= skyFit
		hdr['SKYVAL'] = float(self.sky0)
		# XXX should write fit parameters to header
		if self.skyFits is not None:
			self.skyFits.write(skyFit,extname=extName,header=hdr)
		return data,hdr
	def _finish(self):
		if self.skyFits is not None:
			self.skyFits.close()
			self.skyFits = None

###############################################################################
#                                                                             #
#                               GAIN BALANCE                                  #
#            balances the amplifiers/CCDs with a gain correction              #
#                                                                             #
###############################################################################

class BokCalcGainBalanceFactors(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('read_only',True)
		super(BokCalcGainBalanceFactors,self).__init__(**kwargs)
		self.ampCorStatReg = bokutil.stats_region(kwargs.get('stats_region',
		                                              'amp_corner_ccdcenter'))
		self.statsMethod = kwargs.get('stats_method','mode')
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.clipArgs.setdefault('clip_iters',4)
		self.clipArgs.setdefault('clip_sig',2.5)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.median)
		self.saveArrays = kwargs.get('save_arrays',False)
		self.reset()
	def reset(self):
		self.files = []
		self.gainCors = []
		self.allSkyVals = []
		if self.saveArrays:
			self.arrays = []
	def _preprocess(self,fits,f):
		print 'calculating gain balance factors for ',self.inputNameMap(f)
		self.files.append(f)
		self.skyVals = []
	def _postprocess(self,fits,f):
		skyVals = np.array(self.skyVals)
		# the mean sky level within each CCD
		meanSkyCcd = skyVals.reshape(4,4).mean(axis=-1)
		# the correction needed to balance each amp to the per-CCD sky level
		ampGainCor = np.repeat(meanSkyCcd,4) / skyVals
		# the mean sky level across all CCDs
		meanSky = np.mean(skyVals)
		# the correction needed to balance the CCDs
		ccdGainCor = meanSky / meanSkyCcd
		ccdGainCor = np.repeat(ccdGainCor,4)
		self.gainCors.append((ampGainCor,ccdGainCor))
		self.allSkyVals.append(skyVals)
	def process_hdu(self,extName,data,hdr):
		# XXX would be more efficient to read a subregion
		stats = bokutil.array_stats(data[self.ampCorStatReg],
		                            method=self.statsMethod,
		                            retArray=self.saveArrays,
		                            **self.clipArgs)
		if self.saveArrays:
			skyVal,skyArr = stats
			self.arrays.append(skyArr)
		else:
			skyVal = stats
		self.skyVals.append(skyVal)
		return data,hdr
	def calc_mean_corrections(self):
		gc = np.array(self.gainCors)
		return sigma_clip(gc,iters=2,sig=2.0,axis=0).mean(axis=0).filled()
	def get_values(self):
		return np.array(self.gainCors),np.array(self.allSkyVals)

###############################################################################
#                                                                             #
#                               COMBINE CCDs                                  #
#                balances the amplifiers with a gain correction               #
#                                                                             #
###############################################################################

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

def combine_ccds(fileList,**kwargs):
	inputFileMap = kwargs.get('input_map',bokutil.IdentityNameMap)
	outputFileMap = kwargs.get('output_map',bokutil.IdentityNameMap)
	gainMap = kwargs.get('gain_map')
	origin = kwargs.get('origin','center')
	clobber = kwargs.get('clobber')
	ignoreExisting = kwargs.get('ignore_existing',True)
	# fitsio doesn't accept file descriptors, but tempfile opens files...
	tmpFile = tempfile.NamedTemporaryFile()
	tmpFileName = tmpFile.name
	tmpFile.close()
	# do the extensions in numerical order, instead of HDU list order
	extns = np.array(['IM%d' % ampNum for ampNum in range(1,17)])
	for f in fileList:
		inputFile = inputFileMap(f)
		outputFile = outputFileMap(f)
		print 'combine: ',inputFile,outputFile
		inFits = fitsio.FITS(inputFile)
		if 'CCDJOIN' in inFits[0].read_header():
			print '%s already combined, skipping' % inputFile
			inFits.close()
			continue
		if outputFile != inputFile:
			if os.path.exists(outputFile):
				if clobber:
					os.unlink(outputFile)
				# XXX should be checking header key here?
				elif ignoreExisting:
					print '%s already exists, skipping' % outputFile
					continue
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
			for j,ext in enumerate(extGroup):
				im = inFits[ext].read() 
				hext = inFits[ext].read_header()
				try:
					# copy in the nominal gain values from per-amp headers
					gainKey = 'GAIN%02dA' % int(ext[2:])
					hdr[gainKey] = hext[gainKey]
				except ValueError:
					pass
				if gainMap is not None:
					ampIdx = ampOrder[4*(ccdNum-1)+j] - 1
					gc1,gc2 = gainMap['corrections'][f][:,ampIdx]
					sky = gainMap['skyvals'][f][ampIdx]
					hdr['SKY%02dB'%int(ext[2:])] = sky
					hdr['GAIN%02dB'%int(ext[2:])] = gc1
					if j==0:
						hdr['CCDGAIN'] = gc2
					im *= gc1 * gc2
				ccdIms.append(im)
			# orient the channel images into a mosaic of CCDs and
			# modify WCS & mosaic keywords
			outIm,hdr = _orient_mosaic(hdr,ccdIms,ccdNum,origin)
			outFits.write(outIm,extname='CCD%d'%ccdNum,header=hdr)
		outFits.close()
		if outputFile == inputFile:
			shutil.move(tmpFileName,inputFile)



###############################################################################
#                                                                             #
#                              PROCESS ROUND 2                                #
#   obj detection & masking, divide by supersky flat, ... (CR rejection?)     #
#                                                                             #
###############################################################################

from astropy.convolution.convolve import convolve
from astropy.convolution.kernels import Gaussian2DKernel
from scipy.ndimage.morphology import binary_dilation,binary_closing
import scipy.ndimage.measurements as meas

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
	def __init__(self,**kwargs):
		self.nBin = kwargs.get('binSize',4)
		kwargs.setdefault('header_key','SKYMSK')
		super(BokGenerateSkyFlatMasks,self).__init__(**kwargs)
		self.hiThresh = kwargs.get('high_thresh',5.0)
		self.loThresh = kwargs.get('low_thresh',5.0)
		self.growThresh = kwargs.get('grow_thresh',1.5)
		self.binGrowSize = kwargs.get('mask_grow_size',3)
		self.nSample = kwargs.get('num_sample',4)
		self.nKnots = kwargs.get('num_spline_knots',2)
		self.splineOrder = kwargs.get('spline_order',2)
		self.statsPix = bokutil.stats_region(kwargs.get('stats_region'))
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.clipArgs.setdefault('clip_iters',3)
		self.clipArgs.setdefault('clip_sig',2.2)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.mean)
		self.growKern = None #np.ones((self.binGrowSize,self.binGrowSize),dtype=bool)
	def _preprocess(self,fits,f):
		print 'generating sky mask for ',f
	def process_hdu(self,extName,data,hdr):
		if (data>hdr['SATUR']).sum() > 10000:
			# if too many pixels are saturated mask the whole damn thing
			hdr['BADSKY'] = 1
		n = self.nSample
		sky,rms = bokutil.array_stats(data[self.statsPix][::n,::n],
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
		mask = binary_dilation(mask,mask=(snr>self.growThresh),iterations=0,
		                       structure=self.growKern)
		# fill in holes in the mask
		mask = binary_closing(mask,iterations=5,structure=self.growKern)
		# bright star mask
		bsmask = False # XXX
		# construct the output array
		maskIm = bokutil.magnify(mask,self.nBin) | bsmask
		return maskIm.astype(np.uint8),hdr

def grow_obj_mask(im,objsIm,**kwargs):
	statsPix = bokutil.stats_region(kwargs.get('stats_region',
	                                           'ccd_central_quadrant'))
	growThresh = kwargs.get('grow_thresh',1.0)
	kernelSize = kwargs.get('kernel_size',1.25)
	# XXX missing badpix mask here
	maskedIm = np.ma.masked_array(im,objsIm>0)
	# determine the sky background level and rms
	skym,skys = bokutil.array_stats(maskedIm[statsPix],rms=True)
	# make a pixel-level SNR image
	snrIm = (im - skym) / skys
	# convolve the SNR image to smooth it and slighly grow object footprints
	snrIm = convolve(snrIm,Gaussian2DKernel(kernelSize))
	snrIm[np.isnan(snrIm)] = np.inf
	# grow object mask until pixels reach a threshold in SNR
	mask = binary_dilation(objsIm>0,mask=(snrIm>growThresh),iterations=0)
	# fill holes in the mask
	mask = binary_closing(mask)
	# fix the corners: pre-illumination-correction images have a gradient
	# at the corners, if it is positive (which it is for the upper two CCDs)
	# then objects in those corners are grown until the whole corner is 
	# filled. This simply reverts the corners to the original mask.
	Y,X = np.indices(im.shape)
	corner = Y > ( 2950 + ((2950-4032.)/(4096.-2800))*(X-4096.) )
	mask[corner] = objsIm[corner]
	return mask

def sextract_pass1(fileList,**kwargs):
	clobber = kwargs.get('clobber',False)
	inputNameMap = kwargs.get('input_map',bokutil.IdentityNameMap)
	catalogFileNameMap = kwargs.get('catalog_map',
	                                bokutil.FileNameMap(newSuffix='.cat1'))
	withPsf = kwargs.get('with_psf',False)
	objMaskFileMap = kwargs.get('object_mask_map',
	                             bokutil.FileNameMap(newSuffix='.obj'))
	#bkgImgFileMap = FileNameMap(newSuffix='.back')
	FNULL = open(os.devnull,'w')
	for f in fileList:
		inputFile = inputNameMap(f)
		catalogFile = catalogFileNameMap(f)
		print 'generating ',inputFile,catalogFile
		if os.path.exists(catalogFile) and not clobber:
			continue
		cmd = ['sex','-c',os.path.join(configdir,'bok_pass1.sex'),
		       '-CATALOG_NAME',catalogFile]
		if objMaskFileMap is not None:
			#cmd.extend(['-CHECKIMAGE_TYPE','SEGMENTATION,MINIBACKGROUND',
			#            '-CHECKIMAGE_NAME',
			#                objMaskFileMap(f)+','+bkgImgFileMap(f)])
			cmd.extend(['-CHECKIMAGE_TYPE','SEGMENTATION',
			            '-CHECKIMAGE_NAME','tmpobj.fits'])
#			            '-CHECKIMAGE_NAME',objMaskFileMap(f)])
		cmd.append(inputFile)
		rv = subprocess.call(cmd,stdout=FNULL,stderr=FNULL)
		fits = fitsio.FITS(inputFile,'rw')
		# have to remake the mask file in order to change its data type
		tmpMaskFits = fitsio.FITS('tmpobj.fits')
		if os.path.exists(objMaskFileMap(f)):
			os.unlink(objMaskFileMap(f))
		maskFits = fitsio.FITS(objMaskFileMap(f),'rw')
		maskFits.write(None,header=tmpMaskFits[0].read_header())
		for ccd in ['CCD%d'%i for i in range(1,5)]:
			mask = grow_obj_mask(fits[ccd][:,:],tmpMaskFits[ccd][:,:])
			maskFits.write(mask.astype(np.uint8),extname=ccd,
			               header=tmpMaskFits[ccd].read_header())
		maskFits.close()
		tmpMaskFits.close()
		os.unlink('tmpobj.fits')

class BokNightSkyFlatStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','ccd_central_quadrant')
		kwargs.setdefault('scale','normalize_mode')
		kwargs.setdefault('nsplit',10)
		kwargs.setdefault('fill_value',1.0)
		super(BokNightSkyFlatStack,self).__init__(**kwargs)
		self.clipArgs = { k:v for k,v in kwargs.items() 
		                     if k.startswith('clip_') }
		self.statsMethod = kwargs.get('stats_method','mode')
		self.clipArgs.setdefault('clip_iters',3)
		self.clipArgs.setdefault('clip_sig',4.0)
		self.clipArgs.setdefault('clip_cenfunc',np.ma.mean)
		self.smoothingLength = kwargs.get('smoothing_length',0.05)
		self.rawStackFile = kwargs.get('raw_stack_file')
		self.rawStackFits = None
		self.normCCD = 'CCD1'
		self.headerKey = 'SKYFL'
	def _preprocess(self,fileList,outFits):
		self.norms = np.zeros(len(fileList),dtype=np.float32)
		for i,f in enumerate(fileList):
			fits = bokutil.BokMefImage(self.inputNameMap(f),
			                           mask_file=self.maskNameMap(f),
			                           read_only=True)
			normpix = fits.get(self.normCCD,self.statsPix)
			meanVal = bokutil.array_stats(normpix,method=self.statsMethod,
			                              **self.clipArgs)
			self.norms[i] = 1/meanVal
			print 'norm for image %s is %f' % \
			           (self.inputNameMap(f),meanVal)
		if self.rawStackFile is not None:
			print 'writing raw stack to ',self.rawStackFile(outFits._filename)
			self.rawStackFits = fitsio.FITS(
			                        self.rawStackFile(outFits._filename),'rw')
			# if we've gotten to here, we already know any existing file 
			# needs to be clobbered
			self.rawStackFits.write(None,header=outFits[0].read_header(),
			                        clobber=True)
	def _rescale(self,imCube,scales=None):
		if scales is not None:
			_scales = scales[np.newaxis,:]
		else:
			_scales = self.norms[np.newaxis,:]
		self.scales = _scales.squeeze()
		return imCube * _scales
	def _postprocess(self,extName,stack,hdr):
		# XXX hardcoded params
		stack = interpolate_masked_pixels(stack,along='rows',method='linear')
		# ignore the input mask and adopt the interpolation mask;
		#   nan values mean no interpolation was possible
		interpMask = np.isnan(stack.data)
		cleanStack = np.ma.masked_array(stack.data,mask=interpMask)
		cleanStack = cleanStack.filled(1.0)
		if self.rawStackFile is not None:
			self.rawStackFits.write(cleanStack,extname=extName,header=hdr)
		cleanStack = spline_filter(cleanStack,self.smoothingLength)
		# renormalize to unity, using the combined interp and input mask
		_stack = np.ma.masked_array(cleanStack,mask=interpMask|stack.mask)
		normpix = sigma_clip(_stack[self.statsPix],iters=2,sig=2.5,
		                     cenfunc=np.ma.mean)
		_stack /= normpix.mean()
		return _stack,hdr
	def _cleanup(self):
		super(BokNightSkyFlatStack,self)._cleanup()
		if self.rawStackFits is not None:
			self.rawStackFits.close()
			self.rawStackFits = None



