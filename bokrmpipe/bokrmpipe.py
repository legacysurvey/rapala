#!/usr/bin/env python

import os
import re
import glob
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio

# XXX
from astrotools.idmstuff import loadpath
loadpath()

import bokutil
from bokoscan import BokOverscanSubtract
import bokproc
import badpixels

import boklog
logs = boklog.load_Bok_logs()

datadir = os.environ['HOME']+'/data/observing/Bok/90Prime/RM/'
rdxdir = 'tmprm/'
caldir = rdxdir+'cals/'

if not os.path.exists(rdxdir):
	os.mkdir(rdxdir)
if not os.path.exists(caldir):
	os.mkdir(caldir)

class RMFileNameMap(bokutil.FileNameMap):
	def __init__(self,newSuffix=None,fromRaw=False):
		self.newSuffix = '' if newSuffix is None else newSuffix
		self.fromRaw = fromRaw
	def __call__(self,fileName):
		fn = fileName+self.newSuffix+'.fits'
		if self.fromRaw:
			return os.path.join(datadir,fn+'.gz')
		else:
			return os.path.join(rdxdir,fn)

class MasterBadPixMask(bokutil.FileNameMap):
	def __init__(self,*args,**kwargs):
		# only load it once
		self.fits = fitsio.FITS(os.path.join(caldir,'badpix_master.fits'))
	def __call__(self,*args,**kwargs):
		return self.fits

class MasterBadPixMask4(bokutil.FileNameMap):
	def __init__(self,*args,**kwargs):
		# only load it once
		self.fits = fitsio.FITS(os.path.join(caldir,'badpix_master_4ccd.fits'))
	def __call__(self,*args,**kwargs):
		return self.fits

class processInPlace(object):
	def __init__(self):
		self.fmap = RMFileNameMap()
		self.fremap = {'pass1cat':'.cat1','objmask':'.obj',
		               '_sky':'_tmpsky'}
	def __call__(self,t,output=True):
		if t in self.fremap:
			# some files nned to be remapped even when processing in-place
			return RMFileNameMap(self.fremap[t])
		if output:
			return None
		else:
			return self.fmap

class processToNewFiles(object):
	def __init__(self):
		self.fmap = {'bias':'_b','proc':'_p','comb':'_c',
		             'pass1cat':'.cat1','objmask':'.obj',
		             '_sky':'_tmpsky','sky':'_s','proc2':'_q'}
	def __call__(self,t,output=True):
		return RMFileNameMap(self.fmap[t])

def get_utds(utds=None):
	if utds is None:
		_utds = sorted(logs.keys())
	elif type(utds) is str:
		_utds = [utds,]
	else:
		_utds = utds
	return _utds

def get_files(logs,utds=None,imType=None,filt=None,imRange=None):
	utds = get_utds(utds)
	files = []
	for utd in utds:
		if imRange is None:
			imRange = (0,len(logs[utd]))
		frameNum = np.arange(len(logs[utd]))
		is_range = (imRange[0] <= frameNum) & (frameNum <= imRange[1])
		if imType is None:
			is_type = True
		else:
			is_type = logs[utd]['imType'] == imType
		if filt is None:
			is_filt = True
		else:
			is_filt = logs[utd]['filter'] == filt
		ii = np.where(is_range & is_type & is_filt)[0]
		if len(ii) > 0:
			files.append(char_add('ut'+utd+'/',logs[utd]['fileName'][ii]))
	if len(files)==0:
		return None
	else:
		return np.concatenate(files)

def get_bias_map(utds=None,filt=None):
	biasMap = {}
	biasFiles = sorted(glob.glob(caldir+'bias_*.fits'))
	bias2utd = ".*bias_(\d+).*"
	biasUtds = np.array([int(re.match(bias2utd,fn).groups()[0])
	                       for fn in biasFiles])
	if utds is None:
		utds = sorted(logs.keys())
	for i,utd in enumerate(utds):
		j = np.argmin(np.abs(int(utd)-biasUtds))
		biasFile = biasFiles[j]
		for f in get_files(logs,utd,filt=filt):
			biasMap[f] = biasFile
	return biasMap

def get_flat_map(utds=None,filt=None):
	flatMap = {}
	flatFiles = sorted(glob.glob(caldir+'flat_*.fits'))
	flat2utdfilt = ".*flat_(\d+)_(\w+)_.*"
	utdfilt = [ re.match(flat2utdfilt,fn).groups() for fn in flatFiles]
	flatUtds = np.array([int(utd) for utd,filt in utdfilt])
	flatFilt = np.array([filt for utd,filt in utdfilt])
	if utds is None:
		utds = sorted(logs.keys())
	if filt is None:
		filts = 'gi'
	else:
		filts = filt
	for filt in filts:
		for i,utd in enumerate(utds):
			files = get_files(logs,utd,filt=filt)
			if files is None:
				continue
			jj = np.where(flatFilt==filt)[0]
			j = np.argmin(np.abs(int(utd)-flatUtds[jj]))
			flatFile = flatFiles[jj[j]]
			for f in files:
				flatMap[f] = flatFile
	return flatMap

def overscan_subtract(utds=None,filt=None,**kwargs):
	oscanSubtract = BokOverscanSubtract(input_map=RMFileNameMap(fromRaw=True),
                                        output_map=RMFileNameMap(),
                                        **kwargs)
	oscanSubtract.process_files(get_files(logs,utds,filt=filt))

def make_2d_biases(utds=None,nSkip=2,reject=None,filt=None,**kwargs):
	biasStack = bokproc.BokBiasStack(input_map=RMFileNameMap(),
	                                 reject=reject,
                                     **kwargs)
	for utd in utds:
		files = get_files(logs,utd,imType='zero',filt=filt)
		if files is None:
			continue
		files = files[nSkip:]
		biasNum = 1
		biasStack.stack(files,caldir+'bias_%s_%d.fits' % (utd,biasNum))

def make_dome_flats(file_map,bias_map,utds=None,filt=None,
                    nSkip=1,reject=None,**kwargs):
	bias2Dsub = bokproc.BokCCDProcess(bias_map,
	                                  input_map=RMFileNameMap(),
	                                  output_map=file_map('bias'),
	                                  header_key='BIAS2D',
	                                  **kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=file_map('bias',False),
	                                     **kwargs)
	if filt is None:
		filts = 'gi'
	else:
		filts = filt
	for utd in utds:
		for filt in filts:
			files = get_files(logs,utd,imType='flat',filt=filt)
			if files is None:
				continue
			files = files[nSkip:]
			bias2Dsub.process_files(files)
			flatNum = 1
			flatStack.stack(files,
			                caldir+'flat_%s_%s_%d.fits' % (utd,filt,flatNum))

def make_bad_pixel_masks(**kwargs):
	utd,filt,flatNum = '20140425','g',1
	flatFn = caldir+'flat_%s_%s_%d.fits' % (utd,filt,flatNum)
	bpMaskFile = os.path.join(caldir,'badpix_master.fits')
	badpixels.build_mask_from_flat(flatFn,bpMaskFile,
	             normed_flat_file=flatFn.replace('.fits','_normed.fits'),
	             normed_flat_fit_file=flatFn.replace('.fits','_fit.fits'),
	             binned_flat_file=flatFn.replace('.fits','_binned.fits'),
	             )#,**kwargs)
	bokproc.combine_ccds([bpMaskFile,],
	                     output_map=bokutil.FileNameMap(caldir,'_4ccd'),
	                     apply_gain_correction=False,
	                     **kwargs)

def process_all(file_map,bias_map,flat_map,utds=None,filt=None,
                fixpix=False,**kwargs):
	proc = bokproc.BokCCDProcess(bias_map,
	                             flat_map,
	                             input_map=RMFileNameMap(),
	                             output_map=file_map('proc'),
	                             mask_map=MasterBadPixMask(),
	                             fixpix=fixpix,
	                             **kwargs)
	files = get_files(logs,utds,imType='object',filt=filt)
	proc.process_files(files)
	bokproc.combine_ccds(files,
	                     input_map=file_map('proc',False),
	                     output_map=file_map('comb'),
	                     **kwargs)

def make_supersky_flats(file_map,utds=None,filt=None,
                        skysub=True,**kwargs):
	utds = ['20140427']
	if skysub:
		skySub = bokproc.BokSkySubtract(input_map=file_map('comb',False),
		                                output_map=file_map('_sky'),
		                                mask_map=file_map('objmask'))
		skySub.add_mask(MasterBadPixMask4()())
		stackin = file_map('_sky',False)
	else:
		stackin = file_map('comb',False)
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=file_map('objmask'),
	                    exposure_time_map=bokutil.FileNameMap(caldir,'.exp'),
	                       raw_stack_file=bokutil.FileNameMap(caldir,'_raw'))
	skyFlatStack.set_badpixelmask(MasterBadPixMask4()())
	if filt is None:
		filts = 'gi'
	else:
		filts = filt
	for utd in utds:
		for filt in filts:
			files = get_files(logs,utd,imType='object',filt=filt)
			if files is None:
				continue
			# XXX need to use the bad pix mask for weighting
			bokproc.sextract_pass1(files,
			                       input_map=file_map('comb',False),
			                       catalog_map=file_map('pass1cat'),
			                       object_mask_map=file_map('objmask'),
			                       **kwargs)
			if skysub:
				skySub.process_files(files)
			skyFlatStack.stack(files,
			                   caldir+'skyflat_%s_%s.fits' % (utd,filt))

# process round 2: illum corr and sky sub

def make_images():
	import matplotlib.pyplot as plt
	files = get_files(logs,'20140427',imType='object',filt='g')
	_fmap = RMFileNameMap()
	plt.ioff()
	for ff in files:
		f = _fmap(ff)
		bokproc.make_fov_image_fromfile(f.replace('.fits','_tmpsky.fits'),
		                                f.replace('.fits','.png'),
		                                mask=f.replace('.fits','.obj.fits'))
	plt.ion()

def rmpipe():
	utds = ['20140425','20140427']
	kwargs = {'clobber':False,'verbose':10}
	inplace = True
	filt = 'g'
	# fixpix is sticking nan's into the images in unmasked pixels (???)
	fixpix = False #True
	fileMap = processInPlace() if inplace else processToNewFiles()
	for utd in utds:
		utdir = os.path.join(rdxdir,'ut'+utd)
		if not os.path.exists(utdir): os.mkdir(utdir)
	overscan_subtract(utds,filt=filt,**kwargs)
	make_2d_biases(utds,filt=filt,**kwargs)
	biasMap = get_bias_map(utds,filt=filt)
	make_dome_flats(fileMap,biasMap,utds,filt=filt,**kwargs)
	make_bad_pixel_masks()
	flatMap = get_flat_map(utds,filt=filt)
	process_all(fileMap,biasMap,flatMap,utds,filt=filt,
	            fixpix=fixpix,**kwargs)
	make_supersky_flats(fileMap,utds,filt=filt,**kwargs)
	# XXX for testing
	#fileMap = processToNewFiles()

if __name__=='__main__':
	import sys
	if len(sys.argv)==0:
		rmpipe()
	elif sys.argv[1]=='images':
		make_images()

