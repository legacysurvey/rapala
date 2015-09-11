#!/usr/bin/env python

import os
import fitsio
import numpy as np
from scipy.stats.mstats import mode

from astropy.stats import sigma_clip

def rebin(im,nbin):
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

def stats_region(statreg):
	if statreg is None:
		return np.s_[:,:]
	elif type(statreg) is tuple:
		x1,x2,y1,y2 = statreg
	elif statreg == 'amp_central_quadrant':
		x1,x2,y1,y2 = (512,-512,512,-512)
	elif statreg == 'amp_corner_ccdcenter_small':
		x1,x2,y1,y2 = (-512,-50,-512,-50)
	elif statreg == 'amp_corner_ccdcenter':
		x1,x2,y1,y2 = (-1024,-1,-1024,-1)
	elif statreg == 'centeramp_corner_fovcenter':
		# for the 4 central amps, this is the corner towards the field center
		x1,x2,y1,y2 = (50,1024,50,1024)
	elif statreg == 'ccd_central_quadrant':
		x1,x2,y1,y2 = (1024,-1024,1024,-1024)
	else:
		raise ValueError
	return np.s_[y1:y2,x1:x2]

def _write_stack_header_cards(fileList,cardPrefix):
	hdr = fitsio.read_header(fileList[0])
	for num,f in enumerate(fileList,start=1):
		hdr['%s%03d'%(cardPrefix,num)] = os.path.basename(f)
	return hdr

def build_cube(fileList,extn,masks=None,rows=None):
	if rows is None:
		s = np.s_[:,:]
	else:
		s = np.s_[rows[0]:rows[1],:]
	cube = np.dstack( [ fitsio.FITS(f)[extn][s] for f in fileList ] )
	if masks is not None:
		if isinstance(masks,FileNameMap):
			mask = np.dstack([ fitsio.FITS(masks(f))[extn][s]
			           for f in fileList ])
		else:
			mask = np.dstack([ fitsio.FITS(f)[extn][s] for f in masks ])
	else:
		mask = None
	cube = np.ma.masked_array(cube,mask)
	return cube

class OutputExistsError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class BokMefImage(object):
	'''A wrapper around fitsio that allows the MEF files to be iterated
	   over while updating the data arrays and headers either in-place or
	   to a new file. Also allows for an arbitrary number of masks to be
	   carried with the data.'''
	def __init__(self,fileName,**kwargs):
		self.fileName = fileName
		self.overwrite = kwargs.get('overwrite',False)
		self.keepHeaders = kwargs.get('keep_headers',True)
		self.outFileName = kwargs.get('output_file')
		self.readOnly = kwargs.get('read_only',False)
		self.extensions = kwargs.get('extensions')
		maskFits = kwargs.get('mask_file')
		headerCards = kwargs.get('add_header',{})
		self.closeFiles = []
		if self.readOnly:
			self.fits = fitsio.FITS(fileName)
		elif self.outFileName is None:
			self.fits = fitsio.FITS(fileName,'rw')
			self.outFits = self.fits
			self.clobber = True
			self.closeFiles.append(self.fits)
		else:
			self.fits = fitsio.FITS(fileName)
			if os.path.exists(self.outFileName):
				if self.overwrite:
					os.unlink(self.outFileName)
				else:
					raise OutputExistsError("%s already exists" % 
					                        self.outFileName)
			self.outFits = fitsio.FITS(self.outFileName,'rw')
			hdr = None if not self.keepHeaders else self.fits[0].read_header()
			for k,v in headerCards.items():
				hdr[k] = v
			self.outFits.write(None,header=hdr)
			self.clobber = False
			self.closeFiles.extend([self.fits,self.outFits])
		self.masks = []
		if maskFits is not None:
			self.add_mask(maskFits)
		if self.extensions is None:
			self.extensions = [ h.get_extname().upper() 
			                     for h in self.fits[1:] ]
		self.curExtName = None
	def add_mask(self,maskFits):
		if type(maskFits) is str:
			maskFits = fitsio.FITS(maskFits)
			self.closeFiles.append(maskFits)
		elif type(maskFits) is not fitsio.fitslib.FITS:
			return ValueError
		self.masks.append(maskFits)
	def update(self,data,header=None,noconvert=False):
		if not noconvert:
			# should probably instead track down all the upcasts
			data = data.astype(np.float32)
		self.outFits.write(data,extname=self.curExtName,header=header,
		                   clobber=self.clobber)
	def __iter__(self):
		for self.curExtName in self.extensions:
			data = self.fits[self.curExtName].read()
			hdr = self.fits[self.curExtName].read_header()
			if len(self.masks) > 0:
				mask = self.masks[0][self.curExtName].read().astype(np.bool)
				for m in self.masks[1:]:
					mask |= m[self.curExtName].read().astype(np.bool)
				data = np.ma.masked_array(data,mask=mask)
			yield self.curExtName,data,hdr
	def get(self,extName,subset=None,header=False):
		if subset is None:
			subset = np.s_[:,:]
		data = self.fits[extName][subset]
		if header:
			return data,self.fits[extName].read_header()
		else:
			return data
	def get_header(self,extName):
		return self.fits[extName].read_header()
	def get_xy(self,extName,coordsys='image'):
		hdr = self.fits[extName].read_header()
		return bok_getxy(hdr,coordsys)
	def make_fov_image(self,nbin=1,coordsys='sky'):
		rv = {'coordsys':coordsys,'nbin':nbin}
		hdr0 = self.fits[0].read_header()
		try:
			rv['objname'] = hdr0['OBJECT'].strip()
		except:
			rv['objname'] = 'none'
		for extName,im,hdr in self:
			x,y = bok_getxy(hdr,coordsys)
			if nbin > 1:
				im = rebin(im,nbin)
				x = x[nbin//2::nbin,nbin//2::nbin]
				y = y[nbin//2::nbin,nbin//2::nbin]
			rv[extName] = {'x':x,'y':y,'im':im}
		return rv
	def close(self):
		for fits in self.closeFiles:
			fits.close()

IdentityNameMap = lambda f: f
NullNameMap = lambda f: None

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

class BokProcess(object):
	def __init__(self,**kwargs):
		self.inputNameMap = kwargs.get('input_map',IdentityNameMap)
		self.outputNameMap = kwargs.get('output_map',NullNameMap)
		self.maskNameMap = kwargs.get('mask_map',NullNameMap)
		self.overwrite = kwargs.get('overwrite',False)
		self.keepHeaders = kwargs.get('keep_headers',True)
	def _preprocess(self,fits):
		pass
	def process_hdu(self,extName,data,hdr):
		raise NotImplementedError
	def _postprocess(self):
		pass
	def process_files(self,fileList):
		for f in fileList:
			fits = BokMefImage(self.inputNameMap(f),
			                   output_file=self.outputNameMap(f),
			                   mask_file=self.maskNameMap(f),
			                   keep_headers=self.keepHeaders,
			                   overwrite=self.overwrite)
			self._preprocess(fits)
			for extName,data,hdr in fits:
				data,hdr = self.process_hdu(extName,data,hdr)
				fits.update(data,hdr)
			fits.close()
		self._postprocess()

class BokImArith(BokProcess):
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

class BokMefImageCube(object):
	def __init__(self,**kwargs):
		self.withVariance = kwargs.get('with_variance',False)
		self.scale = kwargs.get('scale')
		self.reject = kwargs.get('reject','sigma_clip')
		self.inputNameMap = kwargs.get('input_map',IdentityNameMap)
		self.maskNameMap = kwargs.get('mask_map',NullNameMap)
		self.statsRegion = kwargs.get('stats_region')
		self.statsPix = stats_region(self.statsRegion)
		self.clipArgs = {'iters':kwargs.get('clip_iters',2),
		                 'sig':kwargs.get('clip_sig',2.5),
		                 'cenfunc':np.ma.mean}
		self.nSplit = kwargs.get('nsplit',1)
		self.overwrite = kwargs.get('overwrite',False)
		self.headerKey = 'CUBE'
		self.extensions = None
	def _rescale(self,imCube,scales=None):
		if scales is not None:
			pass
		elif self.scale is None:
			return imCube
		elif self.scale.startswith('normalize'):
			imScales = imCube[self.statsPix] / imCube[self.statsPix+([0],)]
			imScales = imScales.reshape(-1,imCube.shape[-1])
			scales = sigma_clip(imScales,cenfunc=np.ma.mean,axis=0)
			if self.scale.endswith('_mean'):
				scales = scales.mean(axis=0)
			else:
				# default is the mode
				scales,_ = mode(scales,axis=0)
			scales /= scales.max()
			scales **= -1
		else:
			scales = self.scale(imCube)
		self.scales = scales.squeeze()
		return imCube * scales
	def _reject_pixels(self,imCube):
		if self.reject == 'sigma_clip':
			imCube = sigma_clip(imCube,axis=-1,**self.clipArgs)
		elif self.reject == 'minmax':
			imCube = np.ma.masked_array(imCube)
			imCube[:,:,imCube.argmax(axis=-1)] = np.ma.masked
			imCube[:,:,imCube.argmin(axis=-1)] = np.ma.masked
		return imCube
	def _stack_cube(self,imCube,**kwargs):
		raise NotImplementedError
	def _preprocess(self,fileList):
		pass
	def _postprocess(self,stack,hdr):
		return stack,hdr
	def stack(self,fileList,outputFile,scales=None,**kwargs):
		if os.path.exists(outputFile):
			if self.overwrite:
				os.unlink(outputFile)
			else:
				raise OutputExistsError("%s already exists" % outputFile)
		inputFiles = [ self.inputNameMap(f) for f in fileList ]
		outFits = fitsio.FITS(outputFile,'rw')
		hdr = _write_stack_header_cards(inputFiles,self.headerKey)
		outFits.write(None,header=hdr)
		if self.withVariance:
			varFn = outputFile.replace('.fits','_var.fits')
			try:
				os.unlink(varFn)
			except:
				pass
			varFits = fitsio.FITS(varFn,'rw')
			varFits.write(None,header=hdr)
		extensions = self.extensions
		if extensions is None:
			_fits = fitsio.FITS(inputFiles[0])
			extensions = [ h.get_extname() for h in _fits[1:] ]
		if self.nSplit == 1:
			rowChunks = [None,]
		else:
			# hacky way to divide the array, hardcoded number of rows
			shape = (4032,4096)
			rowSplit = np.arange(0,shape[0],shape[0]//self.nSplit)
			rowSplit[-1] = -1 # grow last split to end of array
			rowChunks = [ (row1,row2) 
			         for row1,row2 in zip(rowSplit[:-1],rowSplit[1:]) ]
		self._preprocess(fileList)
		if self.maskNameMap == NullNameMap:
			# argh, this is a hacky way to check for masks
			masks = None
		else:
			masks = [ self.maskNameMap(f) for f in fileList ]
		for extn in extensions:
			stack = []
			for rows in rowChunks:
				print '::: %s extn %s <%s>' % (outputFile,extn,rows)
				imCube = build_cube(inputFiles,extn,masks=masks,rows=rows)
				imCube = self._rescale(imCube,scales=scales)
				imCube = self._reject_pixels(imCube)
				_stack = self._stack_cube(imCube,**kwargs)
				stack.append(_stack)
			stack = np.vstack(stack)
			hdr = fitsio.read_header(inputFiles[0],extn)
			stack,hdr = self._postprocess(stack,hdr)
			outFits.write(stack.astype(np.float32),extname=extn,header=hdr)
			if self.withVariance:
				var = np.ma.var(imCube,axis=-1).filled(0).astype(np.float32)
				varFits.write(var,extname=extn,header=hdr)
		outFits.close()
		if self.withVariance:
			varFits.close()

class ClippedMeanStack(BokMefImageCube):
	def _stack_cube(self,imCube,weights=None):
		return np.ma.average(imCube,weights=weights,axis=-1)

class MedianStack(BokMefImageCube):
	def _stack_cube(self,imCube,weights=None):
		return np.ma.median(imCube,axis=-1)

