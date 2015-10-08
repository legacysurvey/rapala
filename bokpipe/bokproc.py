#!/usr/bin/env python

import os
import re
import subprocess
import numpy as np
from scipy.stats.mstats import mode
from scipy.interpolate import LSQUnivariateSpline,griddata
from scipy.signal import spline_filter
from astropy.stats import sigma_clip
from astropy.modeling import models,fitting
import fitsio

import bokutil

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

def make_fov_image(fov,pngfn=None,**kwargs):
	import matplotlib.pyplot as plt
	from matplotlib import colors
	from matplotlib import rc
	maskFile = kwargs.get('mask')
	losig = kwargs.get('lo',2.5)
	hisig = kwargs.get('hi',5.0)
	#kwargs.setdefault('cmap',plt.cm.hot_r)
	cmap = plt.cm.jet
	cmap.set_bad('w',1.0)
	w = 0.4575
	h = 0.455
	if maskFile is not None:
		maskFits = fitsio.FITS(maskFile)
	input_vmin = kwargs.get('vmin')
	input_vmax = kwargs.get('vmax')
	rc('text',usetex=False)
	fig = plt.figure(figsize=(6,6.5))
	cax = fig.add_axes([0.1,0.04,0.8,0.01])
	for n,ccd in enumerate(['CCD2','CCD4','CCD1','CCD3']):
		im = fov[ccd]['im']
		if maskFile is not None:
			im = np.ma.masked_array(im,maskFits[ccd][:,:].astype(bool))
		if n == 0:
			i1,i2 = 100//fov['nbin'],1500//fov['nbin']
			if input_vmin is None and input_vmax is None:
				background = sigma_clip(im[i1:i2,i1:i2],iters=3,sig=2.2)
				m,s = background.mean(),background.std()
				print m,s,m-losig*s,m+hisig*s
				vmin = input_vmin if input_vmin is not None else m-losig*s
				vmax = input_vmax if input_vmax is not None else m+hisig*s
			else:
				vmin = input_vmin
				vmax = input_vmax
			norm = colors.Normalize(vmin=vmin,vmax=vmax)
		if im.ndim == 3:
			im = im.mean(axis=-1)
		x = fov[ccd]['x']
		y = fov[ccd]['y']
		i = n % 2
		j = n // 2
		pos = [ 0.0225 + i*w + i*0.04, 0.05 + j*h + j*0.005, w, h ]
		ax = fig.add_axes(pos)
		_im = ax.imshow(im,origin='lower',
		                extent=[x[0,0],x[0,-1],y[0,0],y[-1,0]],
		                norm=norm,cmap=cmap,
		                interpolation=kwargs.get('interpolation','nearest'))
		if fov['coordsys']=='sky':
			ax.set_xlim(x.max(),x.min())
		else:
			ax.set_xlim(x.min(),x.max())
		ax.set_ylim(y.min(),y.max())
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
		if n == 0:
			cb = fig.colorbar(_im,cax,orientation='horizontal')
			cb.ax.tick_params(labelsize=9)
	title = kwargs.get('title',fov.get('file','')+' '+fov.get('objname',''))
	fig.text(0.5,0.99,title,ha='center',va='top',size=12)
	if pngfn is not None:
		plt.savefig(pngfn)
		plt.close(fig)

def make_fov_image_fromfile(fileName,pngfn,nbin=1,coordsys='sky',**kwargs):
	fits = bokutil.BokMefImage(fileName,mask_file=kwargs.get('mask'),read_only=True)
	fov = fits.make_fov_image(nbin,coordsys)
	fov['file'] = fileName
	return make_fov_image(fov,pngfn,**kwargs)

def bok_polyfit(fits,nbin,order,writeImg=False):
	binnedIm = fits.make_fov_image(nbin,'sky')
	# collect the CCD mosaic into a single image
	X,Y,fovIm = [],[],[]
	for ccd in ['CCD%d'%i for i in range(1,5)]:
		print 'getting sky for ',ccd
		clippedIm = sigma_clip(binnedIm[ccd]['im'],iters=2,sig=2.5,
		                       cenfunc=np.ma.mean)
		im = clippedIm.mean(axis=-1)
		#im = binnedIm[ccd]['im'].mean(axis=-1)
		nbad = clippedIm.mask.sum(axis=-1)
		#too_few_pixels = nbad < nbin*2//3
		too_few_pixels = nbad > nbin**2//3
		#print 'biggest: ',nbad.max(),binnedIm[ccd]['im'].shape
		#too_few_pixels = nbad > 1e10
		ii = np.where(~too_few_pixels)
		im[too_few_pixels] = np.ma.masked
		X.append(binnedIm[ccd]['x'][ii])
		Y.append(binnedIm[ccd]['y'][ii])
		fovIm.append(im[ii])
		binnedIm[ccd]['im'] = im
	X = np.concatenate(X)
	Y = np.concatenate(Y)
	fovIm = np.concatenate(fovIm)
	print X.shape,Y.shape,fovIm.shape
	# fit a polynomial to the binned mosaic image
	poly_model = models.Polynomial2D(degree=order)
	#fitfun = fitting.LevMarLSQFitter()
	fitfun = fitting.LinearLSQFitter()
	p = fitfun(poly_model,X,Y,fovIm)
	# return the model images for each CCD at native resolution
	rv = {}
	for ccd in ['CCD%d'%i for i in range(1,5)]:
		rv[ccd] = p(*fits.get_xy(ccd,'sky'))
	rv['skymodel'] = p
	if False: #writeImg:
		# save the original binned image
		make_fov_image(binnedIm,'tmp1.png')
		# and the sky model fit
		for ccd in ['CCD%d'%i for i in range(1,5)]:
			binnedIm[ccd]['im'] = p(binnedIm[ccd]['x'],binnedIm[ccd]['y'])
		make_fov_image(binnedIm,'tmp2.png')
	return rv

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
		flatNorm = mode(stack[self.statsPix],axis=None)[0][0]
		stack /= flatNorm
		try:
			stack = stack.filled(1.0)
		except:
			pass
		hdr['FLATNORM'] = flatNorm
		if self.scales is not None:
			for _i,_scl in enumerate(self.scales,start=1):
				hdr['FLTSCL%02d'%_i] = _scl
		return stack,hdr

class BokNightSkyFlatStack(bokutil.ClippedMeanStack):
	def __init__(self,**kwargs):
		kwargs.setdefault('stats_region','ccd_central_quadrant')
		kwargs.setdefault('scale','normalize_mode')
		kwargs.setdefault('nsplit',10)
		kwargs.setdefault('fill_value',1.0)
		super(BokNightSkyFlatStack,self).__init__(**kwargs)
		self.smoothingLength = kwargs.get('smoothing_length',0.05)
		self.rawStackFile = kwargs.get('raw_stack_file')
		self.rawStackFits = None
		self.normCCD = 'CCD1'
		self.headerKey = 'SKYFLT'
	def _preprocess(self,fileList,outFits):
		self.norms = np.zeros(len(fileList),dtype=np.float32)
		for i,f in enumerate(fileList):
			fits = bokutil.BokMefImage(self.inputNameMap(f),
			                           mask_file=self.maskNameMap(f),
			                           read_only=True)
			normpix = fits.get(self.normCCD,self.statsPix)
			normpix = sigma_clip(normpix,iters=4,sig=2.5,cenfunc=np.ma.mean)
			self.norms[i] = 1/normpix.mean()
			print 'norm for image %s is %f' % \
			           (self.inputNameMap(f),normpix.mean())
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

class BokCCDProcess(bokutil.BokProcess):
	def __init__(self,bias=None,flat=None,**kwargs):
		kwargs.setdefault('header_key','CCDPROC')
		super(BokCCDProcess,self).__init__(**kwargs)
		self.fixPix = kwargs.get('fixpix',False)
		self.fixPixAlong = kwargs.get('fixpix_along','rows')
		self.fixPixMethod = kwargs.get('fixpix_method','linear')
		self.gainMultiply = kwargs.get('gain_multiply',True)
		self.inputGain = kwargs.get('input_gain',{ 'IM%d'%ampNum:g 
		                               for ampNum,g in zip(ampOrder,nominal_gain)})
		#
		self.biasFile = '<none>'
		self.biasFits = None
		self.biasIsMaster = True
		if type(bias) is fitsio.fitslib.FITS:
			self.biasFile = fits._filename
			self.biasFits = bias
		elif type(bias) is str:
			self.biasFile = bias
			self.biasFits = fitsio.FITS(self.biasFile)
		elif type(bias) is dict:
			self.biasMap = bias
			self.biasIsMaster = False
		#
		self.flatFile = '<none>'
		self.flatFits = None
		self.flatIsMaster = True
		if type(flat) is fitsio.fitslib.FITS:
			self.flatFile = fits._filename
			self.flatFits = flat
		elif type(flat) is str:
			self.flatFile = flat
			self.flatFits = fitsio.FITS(self.flatFile)
		elif type(flat) is dict:
			self.flatMap = flat
			self.flatIsMaster = False
	def _preprocess(self,fits,f):
		print 'ccdproc ',fits.fileName,fits.outFileName
		if not self.biasIsMaster:
			biasFile = self.biasMap[f]
			if self.biasFile != biasFile:
				if self.biasFits is not None:
					self.biasFits.close()
				self.biasFile = biasFile
				self.biasFits = fitsio.FITS(self.biasFile)
		if not self.flatIsMaster:
			flatFile = self.flatMap[f]
			if self.flatFile != flatFile:
				if self.flatFits is not None:
					self.flatFits.close()
				self.flatFile = flatFile
				self.flatFits = fitsio.FITS(self.flatFile)
	def process_hdu(self,extName,data,hdr):
		if self.biasFits is not None:
			data -= self.biasFits[extName][:,:]
		if self.flatFits is not None:
			data /= self.flatFits[extName][:,:]
		if self.fixPix:
			data = interpolate_masked_pixels(data,along=self.fixPixAlong,
			                                 method=self.fixPixMethod)
			# what should the fill value be?
			data = data.filled(0)
		hdr['BIASFILE'] = str(self.biasFile)
		hdr['FLATFILE'] = str(self.flatFile)
		if self.fixPix:
			hdr['FIXPIX'] = self.fixPixAlong
		if self.gainMultiply:
			data *= self.inputGain[extName]
			chNum = int(extName.replace('IM',''))
			hdr['GAIN%02dA'%chNum] = self.inputGain[extName]
			hdr['GAIN'] = 1.0 # now in e-
			hdr['SATUR'] = saturation_dn * self.inputGain[extName]
		else:
			hdr['SATUR'] = saturation_dn
		return data,hdr

class BokSkySubtract(bokutil.BokProcess):
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','SKYSUB')
		super(BokSkySubtract,self).__init__(**kwargs)
	def fit_sky_model(self,fits):
		self.skyFit = bok_polyfit(fits,64,1)
		# subtract the sky level at the origin so as to only remove 
		# the gradient
		self.sky0 = self.skyFit['skymodel'](0,0)
	def _preprocess(self,fits,f):
		print 'sky subtracting ',fits.fileName,fits.outFileName
		self.fit_sky_model(fits)
	def process_hdu(self,extName,data,hdr):
		skyFit = (self.skyFit[extName] - self.sky0).astype(np.float32)
		data -= skyFit
		hdr['SKYVAL'] = self.sky0
		skyModel = self.skyFit['skymodel']
		for k,v in zip(skyModel.param_names,skyModel.parameters):
			hdr['SKY'+k.upper()] = v
		return data,hdr

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
		self.skyEst = kwargs.get('sky_estimator',np.ma.mean)
		self.ampCorStatReg = bokutil.stats_region(kwargs.get('stats_region',
		                                                  'amp_corner_ccdcenter'))
		self.clipArgs = {'iters':kwargs.get('clip_iters',3),
		                 'sig':kwargs.get('clip_sig',2.5),
		                 'cenfunc':np.ma.mean}
		self.reset()
	def reset(self):
		self.files = []
		self.gainCors = []
		self.allSkyVals = []
	def _preprocess(self,fits,f):
		print 'calculating gain balance factors for ',f
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
		skyVal = self.skyEst(sigma_clip(data[self.ampCorStatReg],
		                                **self.clipArgs))
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
	inputFileMap = kwargs.get('input_map')
	if inputFileMap is None:
		inputFileMap = bokutil.IdentityNameMap
	outputFileMap = kwargs.get('output_map')
	gainMap = kwargs.get('gain_map')
	origin = kwargs.get('origin','center')
	clobber = kwargs.get('clobber')
	ignoreExisting = kwargs.get('ignore_existing',True)
	tmpFileName = 'tmp.fits'
	# do the extensions in numerical order, instead of HDU list order
	extns = np.array(['IM%d' % ampNum for ampNum in range(1,17)])
	for f in fileList:
		inputFile = inputFileMap(f)
		print 'combine: ',inputFile
		inFits = fitsio.FITS(inputFile)
		if 'CCDJOIN' in inFits[0].read_header():
			print '%s already combined, skipping' % inputFile
			inFits.close()
			continue
		if outputFileMap is not None:
			outFn = outputFileMap(f)
			if os.path.exists(outFn):
				if clobber:
					os.unlink(outFn)
				elif ignoreExisting:
					print '%s already exists, skipping' % outFn
					continue
				else:
					raise bokutil.OutputExistsError('%s already exists'%outFn)
			outFits = fitsio.FITS(outFn,'rw')
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
				# copy in the nominal gain values from per-amp headers
				hext = inFits[ext].read_header()
				gainKey = 'GAIN%02dA' % int(ext[2:])
				hdr[gainKey] = hext[gainKey]
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
		if outputFileMap is None:
			os.rename(tmpFileName,inputFile)



###############################################################################
#                                                                             #
#                              PROCESS ROUND 2                                #
#   obj detection & masking, divide by supersky flat, ... (CR rejection?)     #
#                                                                             #
###############################################################################

from astropy.convolution.convolve import convolve
from astropy.convolution.kernels import Gaussian2DKernel
from scipy.ndimage.morphology import binary_dilation,binary_closing

def grow_obj_mask(im,objsIm,thresh=1.25,**kwargs):
	statsPix = bokutil.stats_region(kwargs.get('stats_region',
	                                           'ccd_central_quadrant'))
	skypix = sigma_clip(im[statsPix],iters=5,sig=2.5,cenfunc=np.ma.mean)
	skym,skys = skypix.mean(),skypix.std()
	snrIm = (im - skym) / skys
	snrIm = convolve(snrIm,Gaussian2DKernel(0.75))
	snrIm[np.isnan(snrIm)] = np.inf
	mask = binary_dilation(objsIm>0,mask=(snrIm>thresh),iterations=0)
	# fill holes on the object footprints
	#mask = binary_closing(mask)
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
	inputNameMap = kwargs.get('input_map')
	if inputNameMap is None:
		inputNameMap = bokutil.IdentityNameMap
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
			#mask = grow_mask(maskFits[ccd][:,:]>0,3)
			mask = grow_obj_mask(fits[ccd][:,:],tmpMaskFits[ccd][:,:])
			#maskFits[ccd].write(mask.astype(np.uint8),clobber=True)
			maskFits.write(mask.astype(np.uint8),extname=ccd,
			               header=tmpMaskFits[ccd].read_header())
		maskFits.close()
		tmpMaskFits.close()
		os.unlink('tmpobj.fits')

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



