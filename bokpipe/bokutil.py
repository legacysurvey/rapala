#!/usr/bin/env python

import os,sys
import types
from time import time
from datetime import datetime
from collections import OrderedDict
import multiprocessing
import fitsio
import numpy as np
from scipy.ndimage.morphology import binary_closing
from astropy.stats import sigma_clip

from bokio import *

# just translates the kwargs
def array_clip(arr,axis=None,**kwargs):
	# for some reason in newer version of astropy (>1.1) axis=-1 
	# no longer works ...
	if axis is not None and axis < 0:
		axis = len(arr.shape) + axis
	arr = sigma_clip(arr,axis=axis,
	                 sigma=kwargs.get('clip_sig',2.5),
	                 iters=kwargs.get('clip_iters',2),
	                 cenfunc=kwargs.get('clip_cenfunc',np.ma.mean))
	return arr

def array_stats(arr,axis=None,method='median',clip=True,rms=False,
                retArray=False,**kwargs):
	if clip:
		arr = array_clip(arr,axis=axis,**kwargs)
	if method=='median':
		val = np.ma.median(arr,axis=axis)
	elif method=='mean':
		val = np.ma.mean(arr,axis=axis)
	elif method=='mode':
		val = 3*np.ma.median(arr,axis=axis) - 2*np.ma.mean(arr,axis=axis)
	else:
		raise ValueError('array stats method %s unrecognized' % method)
	if axis==None:
		# median returns an array but mean doesn't???
		try:
			val = val.compressed()
		except:
			pass
		val = val.squeeze()
	rv = [val]
	if rms:
		_rms = np.ma.std(arr,axis=axis)
		rv.append(_rms)
	if retArray:
		rv.append(arr)
	if len(rv)>1:
		rv = tuple(rv)
	else:
		rv = rv[0]
	return rv

def rebin(im,nbin):
	s = np.array(im.shape) / nbin
	return im.reshape(s[0],nbin,s[1],nbin).swapaxes(1,2).reshape(s[0],s[1],-1)

def magnify(im,nmag):
	n1,n2 = im.shape
	return np.tile(im.reshape(n1,1,n2,1),
	               (1,nmag,1,nmag)).reshape(n1*nmag,n2*nmag)

def bok_getxy(hdr,coordsys='image',coord=None):
	if coord is None:
		y,x = np.indices((hdr['NAXIS2'],hdr['NAXIS1']))
	else:
		x,y = coord
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

# be careful - these slices can't be applied to read a subregion using fitsio.
# i.e., fitsio.FITS(f)[extn][slice] will not work (reads to end of image),
# but fitsio.FITS(f)[extn].read()[slice] will
def stats_region(statreg,stride=None):
	if statreg is None:
		return np.s_[::stride,::stride]
	elif type(statreg) is tuple:
		x1,x2,y1,y2 = statreg
	elif statreg == 'amp_central_quadrant':
		x1,x2,y1,y2 = (512,-512,512,-512)
	elif statreg.startswith('amp_corner_ccdcenter'):
		margin = 5
		sfx = statreg.lstrip('amp_corner_ccdcenter')
		if len(sfx)==0:
			npix = 1024
		else:
			npix = int({'small':512}.get(sfx,sfx))
		x1,x2,y1,y2 = (-npix,-margin,-npix,-margin)
	elif statreg == 'centeramp_corner_fovcenter':
		# for the 4 central amps, this is the corner towards the field center
		x1,x2,y1,y2 = (50,1024,50,1024)
	elif statreg == 'ccd_central_quadrant':
		x1,x2,y1,y2 = (1024,-1024,1024,-1024)
	else:
		raise ValueError
	return np.s_[y1:y2:stride,x1:x2:stride]

def _write_stack_header_cards(fileList,cardPrefix):
	hdr = fitsio.read_header(fileList[0])
	for num,f in enumerate(fileList,start=1):
		hdr['%s%03d'%(cardPrefix,num)] = os.path.basename(f)
	hdr['NCOMBINE'] = len(fileList)
	return hdr

def rows2slice(rows):
	if rows is None:
		s = np.s_[:,:]
	else:
		s = np.s_[rows[0]:rows[1],:]
	return s

def ccd_join(ims,ccdNum,origin='center'):
	'''Given a set of 4 amplifier images (in numerical order, i.e., CCD1 is
	       [IM1,IM2,IM3,IM4]), combine them into a single CCD image.
	     origin='center' means (x,y) = (0,0) is at the corner nearest to
	            the center of the focal plane.
	'''
	im1,im2,im3,im4 = ims
	# if the arrays are masked, np.[hv]stack will unmask all elements.
	# but if they are not masked, np.ma.[hv]stack will make them masked.
	if any([isinstance(im,np.ma.masked_array) for im in ims]):
		hstack,vstack = np.ma.hstack,np.ma.vstack
	else:
		hstack,vstack = np.hstack,np.vstack
	# orient the channel images N through E and stack into CCD image
	ccdIm = vstack([ hstack([ np.flipud(np.rot90(im2)),
	                          np.rot90(im4,3) ]),
	                 hstack([ np.rot90(im1),
	                          np.fliplr(np.rot90(im3)) ]) ])
	if origin == 'lower left':
		pass
	elif origin == 'center':
		if ccdNum == 1:
			ccdIm = np.fliplr(ccdIm)
		elif ccdNum == 2:
			ccdIm = np.rot90(ccdIm,2)
		elif ccdNum == 3:
			pass
		elif ccdNum == 4:
			ccdIm = np.flipud(ccdIm)
	return ccdIm

def ccd_split(ccdIm,ccdNum,origin='center'):
	'''Inverse of ccd_join: given a CCD image, split it into 4 amplifier
	   images with their original orientation. Returned in numerical order
	   (i.e., [IM1,IM2,IM3,IM4] for CCD1)'''
	# see above (ccd_join).
	if origin == 'lower left':
		pass
	elif origin == 'center':
		if ccdNum == 1:
			ccdIm = np.fliplr(ccdIm)
		elif ccdNum == 2:
			ccdIm = np.rot90(ccdIm,2)
		elif ccdNum == 3:
			pass
		elif ccdNum == 4:
			ccdIm = np.flipud(ccdIm)
	tmp1,tmp2 = np.vsplit(ccdIm,2)
	im2,im4 = np.hsplit(tmp1,2)
	im1,im3 = np.hsplit(tmp2,2)
	im2 = np.rot90(np.flipud(im2),-1)
	im4 = np.rot90(im4,-3)
	im1 = np.rot90(im1,-1)
	im3 = np.rot90(np.fliplr(im3),-1)
	return [im1,im2,im3,im4]

def build_cube(fileList,extn,masks=None,rows=None,badKey=None,
               maskType='gtzero'):
	s = rows2slice(rows)
	cube = np.dstack( [ fitsio.FITS(f)[extn][s] for f in fileList ] )
	_masks = []
	if masks is not None:
		if isinstance(masks,FileNameMap):
			maskFiles = [ masks(f) for f in fileList ]
		else:
			maskFiles = masks
		for f in maskFiles:
			hdu = fitsio.FITS(f)[extn]
			_masks.append(load_mask(hdu[s],maskType))
			# hacky to put this special case here...
			if badKey is not None:
				hdr = hdu.read_header()
				if badKey in hdr:
					_masks[-1][:] = True
		mask = np.dstack(_masks).astype(np.bool)
	else:
		mask = None
	cube = np.ma.masked_array(cube,mask)
	return cube

class OutputExistsError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

def get_timestamp():
	return datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')

class TimerLog():
	def __init__(self):
		self.stages = ['Start']
		self.times = [time()]
	def __call__(self,stage):
		self.stages.append(stage)
		self.times.append(time())
	def dump(self):
		self.__call__('Finish')
		stages = self.stages[1:]
		times = np.array(self.times[1:]) - self.times[0]
		itimes = np.diff(self.times)
		ftimes = itimes / times[-1]
		tscale = 1.
		tunit = 'sec'
		if times[-1] > 5*3600:
			tscale = 3600.  # -> hours
			tunit = 'hrs'
		elif times[-1] > 5*60:
			tscale = 60.  # -> minutes
			tunit = 'min'
		itimes /= tscale
		times /= tscale
		print '%20s %8s %8s %8s [%s]' % ('stage','time','elapsed','frac',tunit)
		for t in zip(stages,itimes,times,ftimes):
			print '%20s %8.3f %8.3f %8.3f' % t
		print

def mask_saturation(extName,data,correct_inverted=True):
	satVal = {'IM5':55000,'IM7':55000}.get(extName,62000)
	mask = data > satVal
	if correct_inverted:
		filledmask = binary_closing(mask,iterations=20)
		data[filledmask&~mask] = 65535
		mask = filledmask
	return data,mask

def load_mask(maskIm,maskType):
	if maskIm.dtype == np.dtype(np.bool):
		return maskIm
	elif maskType=='gtzero':
		return maskIm > 0
	elif maskType=='nonzero':
		return maskIm != 0
	elif maskType=='zerobad':
		return maskIm == 0
	else:
		raise ValueError

class BokMefImage(object):
	'''A wrapper around fitsio that allows the MEF files to be iterated
	   over while updating the data arrays and headers either in-place or
	   to a new file. Also allows for an arbitrary number of masks to be
	   carried with the data.'''
	def __init__(self,fileName,**kwargs):
		self.fileName = fileName
		self.outFileName = kwargs.get('output_file')
		self.clobber = kwargs.get('clobber',False)
		self.headerKey = kwargs.get('header_key')
		self.readOnly = kwargs.get('read_only',False)
		self.keepHeaders = kwargs.get('keep_headers',True)
		self.extensions = kwargs.get('extensions')
		maskFits = kwargs.get('mask_file')
		maskType = kwargs.get('mask_type','gtzero')
		headerCards = kwargs.get('add_header',{})
		self.headerFixes = kwargs.get('header_fixes',[])
		self.closeFiles = []
		if self.readOnly:
			self.fits = fitsio.FITS(self.fileName)
		else:
			if self.outFileName == self.fileName:
				self._check_header_key(self.fileName)
				self.outFits = self.fits = fitsio.FITS(self.fileName,'rw')
				self.closeFiles.append(self.fits)
				self.clobberHdus = True
				for extNum,hdrfix in self.headerFixes:
					if extNum == 0:
						for k,v in hdrfix:
							self.outFits[0].write_key(k,v)
				if self.headerKey is not None:
					self.outFits[0].write_key(self.headerKey,get_timestamp())
			else:
				if os.path.exists(self.outFileName):
					# first see if the output file has already generated
					if not self.clobber:
						self._check_header_key(self.outFileName)
					# can't seem to overwrite extension 0 with fitsio, so
					# for now just deleting the existing file
					os.unlink(self.outFileName)
				self.clobberHdus = False
				self.fits = fitsio.FITS(self.fileName)
				self.outFits = fitsio.FITS(self.outFileName,'rw')
				self.closeFiles.extend([self.fits,self.outFits])
				if self.keepHeaders:
					hdr = self.fits[0].read_header()
				else:
					hdr = {}
				for k,v in headerCards.items():
					hdr[k] = v
				for extNum,hdrfix in self.headerFixes:
					if extNum == 0:
						for k,v in hdrfix.items():
							hdr[k] = v
				if self.headerKey is not None:
					hdr[self.headerKey] = get_timestamp()
				self.outFits.write(None,header=hdr)
		self.masks = []
		self.maskTypes = []
		if maskFits is not None:
			if not isinstance(maskFits,FakeFITS):
				try:
					m = maskFits(fileName)
				except:
					m = FakeFITS(maskFits)
			self.add_mask(maskFits,maskType)
		if self.extensions is None:
			self.extensions = [ h.get_extname().upper() 
			                     for h in self.fits[1:] ]
		self.curExtName = None
	def _check_header_key(self,fileName):
		if self.headerKey is not None:
			hdr0 = fitsio.read_header(fileName,0)
			if self.headerKey in hdr0:
				raise OutputExistsError('key %s already exists' % 
				                        self.headerKey)
	def add_mask(self,maskFits,maskType='gtzero'):
		if isinstance(maskFits,basestring):
			maskFits = fitsio.FITS(maskFits)
			self.closeFiles.append(maskFits)
		#elif not isinstance(maskFits,fitsio.FITS):
		#	return ValueError
		self.masks.append(maskFits)
		self.maskTypes.append(maskType)
	def update(self,data,header=None,noconvert=False):
		if self.readOnly:
			return
		if not noconvert:
			# should probably instead track down all the upcasts
			data = data.astype(np.float32)
		# apply any header fixes
		if header is not None:
			for extNum,hdrfix in self.headerFixes:
				if extNum == self.curExtName:
					for k,v in hdrfix.items():
						header[k] = v
			# NOAO archive adds these keywords
			header.delete('BZERO')
			header.delete('BSCALE')
		# I thought this was overwriting existing HDUs, but doesn't seem to..
		#self.outFits.write(data,extname=self.curExtName,header=header,
		#                   clobber=self.clobberHdus)
		if self.clobberHdus:
			self.outFits[self.curExtName].write(data)
			self.outFits[self.curExtName].write_keys(header)
		else:
			self.outFits.write(data,extname=self.curExtName,header=header,
			                   clobber=False)
	def _load_masks(self,extName,subset):
		if subset is None:
			subset = np.s_[:,:]
		mask = load_mask(self.masks[0][extName][subset],self.maskTypes[0])
		for m,mtyp in zip(self.masks[1:],self.maskTypes[1:]):
			mask |= load_mask(m[extName][subset],mtyp)
		return mask
	def __iter__(self):
		for self.curExtName in self.extensions:
			data = self.fits[self.curExtName].read()
			hdr = self.fits[self.curExtName].read_header()
			if len(self.masks) > 0:
				mask = self._load_masks(self.curExtName,None)
				data = np.ma.masked_array(data,mask=mask)
			yield self.curExtName,data,hdr
	def get(self,extName,subset=None,header=False):
		if subset is None:
			subset = np.s_[:,:]
		data = self.fits[extName][subset]
		if len(self.masks) > 0:
			mask = self._load_masks(extName,subset)
			data = np.ma.masked_array(data,mask=mask)
		if header:
			return data,self.fits[extName].read_header()
		else:
			return data
	def get_header(self,extName):
		return self.fits[extName].read_header()
	def get_xy(self,extName,coordsys='image'):
		hdr = self.fits[extName].read_header()
		return bok_getxy(hdr,coordsys)
	def make_fov_image(self,nbin=1,coordsys='sky',
	                   binfunc=None,binclip=False,single=False,mingood=0):
		rv = OrderedDict()
		hdr0 = self.fits[0].read_header()
		for extName,im,hdr in self:
			x,y = bok_getxy(hdr,coordsys)
			if nbin > 1:
				im = rebin(im,nbin)
				if binclip:
					im = array_clip(im,axis=-1)
				if mingood > 0:
					nbad = im.mask.sum(axis=-1)
				if binfunc is not None:
					im = binfunc(im,axis=-1)
				if mingood > 0:
					im[nbad>mingood] = np.ma.masked
				x = x[nbin//2::nbin,nbin//2::nbin]
				y = y[nbin//2::nbin,nbin//2::nbin]
			rv[extName] = {'x':x,'y':y,'im':im}
		if single:
			_rv = {}
			for f in ['x','y','im']:
				dstack = np.ma.dstack if f=='im' else np.dstack
				_rv[f] = dstack([rv[extName][f] for extName in rv])
			rv = _rv
		# add the object name and binning parameters
		try:
			rv['objname'] = hdr0['OBJECT'].strip()
		except:
			rv['objname'] = 'none'
		rv['coordsys'] = coordsys
		rv['nbin'] = nbin
		return rv
	def close(self):
		for fits in self.closeFiles:
			fits.close()

# make the instance methods pickleable using code from 
# https://gist.github.com/bnyeggen/1086393

def _pickle_method(method):
	func_name = method.im_func.__name__
	obj = method.im_self
	cls = method.im_class
	if func_name.startswith('__') and not func_name.endswith('__'): #deal with mangled names
		cls_name = cls.__name__.lstrip('_')
		func_name = '_' + cls_name + func_name
	return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
	for cls in cls.__mro__:
		try:
			func = cls.__dict__[func_name]
		except KeyError:
			pass
		else:
			break
	return func.__get__(obj, cls)

import copy_reg
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

def mplog(msg,nProc=None):
	if nProc != 1:
		pname = multiprocessing.current_process().name
		if not (nProc is None and pname == 'MainProcess'):
			pid = pname.split('-')[1]
			print '[%2s] '%pid,
	print msg

class FakeFITS(object):
	'''Make a fake fitsio.FITS class that stores the image data from a FITS
	   object in a dictionary. This object is pickle-able and so can be used
	   with the multiprocessing feature in BokProcess, unlike fitsio.FITS.'''
	def __init__(self,fits):
		if isinstance(fits,basestring):
			fits = fitsio.FITS(fits)
		# this gets accessed for FITS objects, so need to replicate it
		self._filename = fits._filename
		self.data = [None] # empty first extension, like FITS MEF
		self.extMap = {}
		for extNum,hdu in enumerate(fits[1:],start=1):
			self.data.append(hdu.read())
			self.extMap[hdu.get_extname().upper()] = extNum
	def __getitem__(self,extn):
		'''index either by extension number or name'''
		if isinstance(extn,basestring):
			extn = self.extMap[extn.upper()]
		return self.data[extn]
	def close(self):
		'''need to provide hook for this because fits.close() is used often'''
		pass

class BokProcess(object):
	_procMsg = '<BokProcess> %s'
	def __init__(self,**kwargs):
		self.inputNameMap = kwargs.get('input_map',IdentityNameMap)
		self.outputNameMap = kwargs.get('output_map',IdentityNameMap)
		self.masks = []
		self.maskTypes = []
		self.maskNameMap = kwargs.get('mask_map')
		self.maskType = kwargs.get('mask_type','gtzero')
		if self.maskNameMap is None:
			self.maskNameMap = NullNameMap
		elif isinstance(self.maskNameMap,fitsio.FITS):
			# a master mask instead of a map
			self.add_mask(self.maskNameMap,self.maskType)
			self.maskNameMap = NullNameMap
		self.clobber = kwargs.get('clobber',False)
		self.readOnly = kwargs.get('read_only',False)
		self.headerKey = kwargs.get('header_key')
		self.headerFixes = kwargs.get('header_fixes',{})
		self.ignoreExisting = kwargs.get('ignore_existing',True)
		self.keepHeaders = kwargs.get('keep_headers',True)
		self.extensions = kwargs.get('extensions')
		self.verbose = kwargs.get('verbose',0)
		self.debug = kwargs.get('debug',False)
		self.nProc = kwargs.get('processes',1)
		self.procMap = kwargs.get('procmap',map)
		self.noConvert = False
	def add_mask(self,maskFits,maskType='gtzero'):
		if not isinstance(maskFits,FakeFITS):
			try:
				maskFits = maskFits(None) # assume it is a map
			except:
				maskFits = FakeFITS(maskFits)
		self.masks.append(maskFits)
		self.maskTypes.append(maskType)
	def _proclog(self,f):
		mplog(self._procMsg % f, self.nProc)
	def _preprocess(self,fits,f):
		self._proclog(f)
	def process_hdu(self,extName,data,hdr):
		raise NotImplementedError
	def _postprocess(self,fits,f):
		pass
	def _finish(self):
		return None
	def _getOutput(self):
		return None
	def _ingestOutput(self,procOutput):
		'''parallel routines may need to collect output from processing
	       of individual files before finishing.'''
		pass
	def process_file(self,f):
		try:
			fits = BokMefImage(self.inputNameMap(f),
			                   output_file=self.outputNameMap(f),
			                   mask_file=self.maskNameMap(f),
			                   mask_type=self.maskType,
			                   keep_headers=self.keepHeaders,
			                   clobber=self.clobber,
			                   header_key=self.headerKey,
			                   header_fixes=self.headerFixes.get(f,{}),
			                   read_only=self.readOnly,
			                   extensions=self.extensions)
		except OutputExistsError,msg:
			if self.ignoreExisting:
				if self.verbose > 0:
					_f = self.outputNameMap(f)
					print '%s already processed by %s'%(_f,self.headerKey)
				return
			else:
				raise OutputExistsError(msg)
		for maskIm,maskType in zip(self.masks,self.maskTypes):
			fits.add_mask(maskIm,maskType)
		self._preprocess(fits,f)
		for extName,data,hdr in fits:
			data,hdr = self.process_hdu(extName,data,hdr)
			fits.update(data,hdr,noconvert=self.noConvert)
		self._postprocess(fits,f)
		fits.close()
		return self._getOutput()
	def _null_result(self,f):
		return None
	def _process_file_exc(self,f):
		if self.debug:
			return self.process_file(f)
		else:
			try:
				return self.process_file(f)
			except Exception,e:
				sys.stderr.write('ERROR: failed to process %s [%s]\n' %
				                       (self.inputNameMap(f),e))
				return self._null_result(f)
	def _process_file_group(self,fgrp):
		return map(self._process_file_exc,fgrp)
	def process_files(self,fileList):
		# pool objects can't be pickled so have to save it, remove it from
		# object, then restore it
		procMap = self.procMap
		self.procMap = None
		if isinstance(fileList[0],types.StringTypes):
			procOut = procMap(self._process_file_exc,fileList)
		else:
			procOut = procMap(self._process_file_group,fileList)
			procOut = [ out for outgrp in procOut for out in outgrp ]
		procOut = filter(lambda p: p is not None,procOut)
		if self.nProc > 1:
			self._ingestOutput(procOut)
		self.procMap = procMap
		self._finish()

class BokMefImageCube(object):
	def __init__(self,**kwargs):
		self.withVariance = kwargs.get('with_variance',False)
		self.scale = kwargs.get('scale')
		self.reject = kwargs.get('reject','sigma_clip')
		self.inputNameMap = kwargs.get('input_map',IdentityNameMap)
		self.outputNameMap = kwargs.get('output_map',IdentityNameMap)
		self.maskNameMap = kwargs.get('mask_map',NullNameMap)
		self.maskType = kwargs.get('mask_type','gtzero')
		self.badKey = kwargs.get('header_bad_key')
		self.expTimeNameMap = kwargs.get('exposure_time_map',NullNameMap)
		self.withExpTimeMap = self.expTimeNameMap != NullNameMap
		self.statsRegion = kwargs.get('stats_region')
		self.statsStride = kwargs.get('stats_stride')
		self.statsPix = stats_region(self.statsRegion,self.statsStride)
		self.clipArgs = {'iters':kwargs.get('clip_iters',2),
		                 'sig':kwargs.get('clip_sig',2.5),
		                 'cenfunc':kwargs.get('clip_cenfunc',np.ma.mean)}
		self.fillValue = kwargs.get('fill_value',np.nan)
		self.maxMemBytes = self.maxMemGB = kwargs.get('maxmem')
		if self.maxMemBytes:
			self.maxMemBytes *= 1024**3
		self.clobber = kwargs.get('clobber',False)
		self.ignoreExisting = kwargs.get('ignore_existing',True)
		self.deleteFiles = kwargs.get('delete_files',False)
		self.verbose = kwargs.get('verbose',0)
		self.headerKey = 'CUBE'
		self.extensions = None
		self.badPixelMask = None
		self.scaleKey = kwargs.get('scale_key','IMSCL')
		self._scales = None
		self.minNexp = None
	def set_badpixelmask(self,maskFits):
		if isinstance(maskFits,FakeFITS):
			self.badPixelMask = maskFits
		else:
			# hope that either it is a map object or an already-loaded FITS...
			try:
				self.badPixelMask = maskFits(None)
			except:
				self.badPixelMask = FakeFITS(maskFits)
	def _rescale(self,imCube,scales=None):
		if scales is not None:
			pass
		elif self.scale is None:
			return imCube
		elif self._scales is not None:
			scales = self._scales
		elif self.scale.startswith('normalize'):
			method = self.scale[self.scale.find('_')+1:]
			imScales = imCube[self.statsPix] / imCube[self.statsPix+([0],)]
			imScales = imScales.reshape(-1,imCube.shape[-1])
			scales = array_stats(imScales,axis=0,method=method)
			scales /= scales.max()
			scales **= -1
		else:
			scales = self.scale(imCube)
		self.scales = scales.squeeze()
		# save the scales that were used. note that for nsplit>1, this has
		# the effect of using the scales from the first chunk only.
		self._scales = scales
		return imCube * scales
	def _reject_pixels(self,imCube):
		if self.reject == 'sigma_clip':
			imCube = array_clip(imCube,axis=-1,**self.clipArgs)
		elif self.reject == 'minmax':
			imCube = np.ma.masked_array(imCube)
			imCube[:,:,imCube.argmax(axis=-1)] = np.ma.masked
			imCube[:,:,imCube.argmin(axis=-1)] = np.ma.masked
		return imCube
	def _load_weights(self,weights,fileList,extn,rows):
		# if it's a map convert it to a list of files
		if isinstance(weights,FileNameMap):
			weights = [weights(f) for f in fileList]
		# if it's a list of files convert it to arrays
		if type(weights) is list:
			weights = build_cube(weights,extn,rows=rows)
		# return either the arrays, and input weight array, or None
		return weights
	def _stack_cube(self,imCube,weights=None,**kwargs):
		raise NotImplementedError
	def _preprocess(self,fileList,outFits):
		pass
	def _postprocess(self,extName,stack,hdr):
		if self._scales is not None:
			for i,s in enumerate(self._scales):
				hdr[self.scaleKey+'%03d'%i] = float(s)
		return stack,hdr
	def stack(self,fileList,outputFile,weights=None,scales=None,**kwargs):
		outputFile = self.outputNameMap(outputFile)
		if os.path.exists(outputFile):
			if self.clobber:
				clobberHdus = True
			else:
				if self.ignoreExisting:
					if self.verbose > 0:
						print '%s already stacked' % outputFile
					return
				else:
					raise OutputExistsError("%s already exists" % outputFile)
		else:
			clobberHdus = False
		inputFiles = map(self.inputNameMap,fileList)
		outFits = fitsio.FITS(outputFile,'rw',clobber=clobberHdus)
		hdr = _write_stack_header_cards(inputFiles,self.headerKey)
		outFits.write(None,header=hdr)
		if self.withExpTimeMap:
			expFn = self.expTimeNameMap(outputFile)
			try:
				os.unlink(expFn)
			except:
				pass
			expTimeFits = fitsio.FITS(expFn,'rw')
			expTimeFits.write(None,header=hdr)
			expTimes = [ fitsio.read_header(_f,ext=0)['EXPTIME']
			                 for _f in inputFiles]
			expTimes = np.array(expTimes).astype(np.float32)
			expTimes = expTimes[np.newaxis,np.newaxis,:]
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
		if self.maxMemBytes is None:
			nsplits = 1
		else:
			numRows,numCols = 4032,4096
			nbytes = np.dtype(np.float32).itemsize
			stacksize = nbytes * numRows * numCols * len(inputFiles)
			nsplits = stacksize // int(self.maxMemBytes) + 1
		if nsplits == 1:
			rowChunks = [None,]
		else:
			rowsplits = np.linspace(0,numRows,nsplits+1,dtype=np.int32)
			rowChunks = [ (row1,row2) 
			         for row1,row2 in zip(rowsplits[:-1],rowsplits[1:]) ]
		self._preprocess(fileList,outFits)
		if self.maskNameMap == NullNameMap:
			# argh, this is a hacky way to check for masks
			masks = None
		else:
			masks = [ self.maskNameMap(f) for f in fileList ]
		for extn in extensions:
			stack = []
			if self.withExpTimeMap:
				expTime = []
			if self.withVariance:
				var = []
			for rows in rowChunks:
				print '::: %s extn %s <%s>' % (outputFile,extn,rows)
				imCube = build_cube(inputFiles,extn,masks=masks,rows=rows,
				                    badKey=self.badKey,maskType=self.maskType)
				imCube = self._rescale(imCube,scales=scales)
				imCube = self._reject_pixels(imCube)
				w = self._load_weights(weights,fileList,extn,rows)
				_stack = self._stack_cube(imCube,w,**kwargs)
				if self.badPixelMask is not None:
					bpmsk = self.badPixelMask[extn][rows2slice(rows)]
					_stack.mask |= load_mask(bpmsk,'gtzero')
				if self.minNexp is not None:
					nexp = np.sum(~imCube.mask,axis=-1)
					_stack.mask |= nexp < self.minNexp
				stack.append(_stack)
				if self.withExpTimeMap:
					expTime.append(np.sum(~imCube.mask*expTimes,axis=-1))
				if self.withVariance:
					# XXX this isn't the right variance for a weighted sum,
					#     really the var calculation needs to happen in 
					#     _stack_cube since it is implementation-dependent
					var.append(np.ma.var(imCube,axis=-1))
			stack = np.ma.vstack(stack)
			hdr = fitsio.read_header(inputFiles[0],extn)
			stack,hdr = self._postprocess(extn,stack,hdr)
			try:
				finalStack = stack.filled(self.fillValue).astype(np.float32)
			except AttributeError:
				finalStack = stack.astype(np.float32)
			outFits.write(finalStack,extname=extn,header=hdr)
			if self.withExpTimeMap:
				expTime = np.ma.vstack(expTime)
				expTimeFits.write(expTime,extname=extn,header=hdr)
			if self.withVariance:
				var = np.ma.vstack(var)
				var = var.filled(0).astype(np.float32)
				varFits.write(var,extname=extn,header=hdr)
		outFits.close()
		if self.withExpTimeMap:
			expTimeFits.close()
		if self.withVariance:
			varFits.close()
		if self.deleteFiles:
			map(os.unlink,inputFiles)
		self._cleanup()
	def _cleanup(self):
		self._scales = None

class ClippedMeanStack(BokMefImageCube):
	def _stack_cube(self,imCube,weights=None):
		# why does it get promoted?
		return np.ma.average(imCube,weights=weights,axis=-1).astype(np.float32)

class MedianStack(BokMefImageCube):
	def _stack_cube(self,imCube,weights=None):
		return np.ma.median(imCube,axis=-1).astype(np.float32)

