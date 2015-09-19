#!/usr/bin/env python

import os
import re
import glob
import numpy as np
from numpy.core.defchararray import add as char_add

# XXX
from astrotools.idmstuff import loadpath
loadpath()

import bokutil
from bokoscan import BokOverscanSubtract
import bokproc

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

class processInPlace(object):
	def __init__(self):
		self.fmap = RMFileNameMap()
	def __call__(self,t,output=True):
		if output:
			return None
		else:
			return self.fmap

class processToNewFiles(object):
	def __init__(self):
		self.fmap = {'bias':'_b','proc':'_p','comb':'_c',
		             'pass1cat':'.cat1','objmask':'.obj',
		             'sky':'_s','proc2':'_q'}
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

def get_bias_map(utds=None):
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
		for f in get_files(logs,utd):
			biasMap[f] = biasFile
	return biasMap

def get_flat_map(utds=None):
	flatMap = {}
	flatFiles = sorted(glob.glob(caldir+'flat_*.fits'))
	flat2utdfilt = ".*flat_(\d+)_(\w+)_.*"
	utdfilt = [ re.match(flat2utdfilt,fn).groups() for fn in flatFiles]
	flatUtds = np.array([int(utd) for utd,filt in utdfilt])
	flatFilt = np.array([filt for utd,filt in utdfilt])
	if utds is None:
		utds = sorted(logs.keys())
	for filt in 'gi':
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

def overscan_subtract(utds=None,**kwargs):
	oscanSubtract = BokOverscanSubtract(input_map=RMFileNameMap(fromRaw=True),
                                        output_map=RMFileNameMap(),
                                        **kwargs)
	oscanSubtract.process_files(get_files(logs,utds))

def make_2d_biases(utds=None,nSkip=2,reject=None,**kwargs):
	biasStack = bokproc.BokBiasStack(input_map=RMFileNameMap(),
	                                 reject=reject,
                                     **kwargs)
	for utd in utds:
		files = get_files(logs,utd,imType='zero')
		if files is None:
			continue
		files = files[nSkip:]
		biasNum = 1
		biasStack.stack(files,caldir+'bias_%s_%d.fits' % (utd,biasNum))

def make_dome_flats(file_map,bias_map,utds=None,nSkip=1,reject=None,**kwargs):
	bias2Dsub = bokproc.BokDebiasFlatten(bias_map,
	                                     input_map=RMFileNameMap(),
	                                     output_map=file_map('bias'),
	                                     header_key='BIAS2D',
	                                     **kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=file_map('bias',False),
	                                     **kwargs)
	for utd in utds:
		for filt in 'gi':
			files = get_files(logs,utd,imType='flat',filt=filt)
			if files is None:
				continue
			files = files[nSkip:]
			print filt,files
			bias2Dsub.process_files(files)
			flatNum = 1
			flatStack.stack(files,
			                caldir+'flat_%s_%s_%d.fits' % (utd,filt,flatNum))

def process_all(file_map,bias_map,flat_map,utds=None,**kwargs):
	proc = bokproc.BokDebiasFlatten(bias_map,
	                                flat_map,
	                                input_map=RMFileNameMap(),
	                                output_map=file_map('proc'),
	                                **kwargs)
	files = get_files(logs,utds,imType='object')
	proc.process_files(files)
	bokproc.combine_ccds(files,
	                     input_map=file_map('proc',False),
	                     output_map=file_map('comb'),
	                     **kwargs)

def make_supersky_flats(file_map,utds=None,skysub=False,**kwargs):
	utds = ['ut20140427']
	if skysub:
		skySub = bokproc.BokSkySubtract(input_map=file_map('comb',False),
		                                output_map=file_map('sky'),
		                                mask_map=file_map('objmask'))
		stackin = file_map('sky',False)
	else:
		stackin = file_map('comb',False)
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=file_map('objmask'))
	for utd in utds:
		for filt in 'gi':
			files = get_files(logs,utd,imType='object',filt=filt)
			bokproc.sextract_pass1(files,
			                       input_map=file_map('comb',False),
			                       catalog_mask_map=file_map('pass1cat'),
			                       object_mask_map=file_map('objmask'),
			                       **kwargs)
			if skysub:
				skySub.process_files(files)
			skyFlatStack.stack(files,
			                   caldir+'skyflat_%s_%s.fits' % (utd,filt))

# process round 2: illum corr and sky sub

def rmpipe():
	utds = ['20140425','20140427']
	kwargs = {'clobber':False,'verbose':10}
	inplace = True
	fileMap = processInPlace() if inplace else processToNewFiles()
	for utd in utds:
		utdir = os.path.join(rdxdir,'ut'+utd)
		if not os.path.exists(utdir): os.mkdir(utdir)
	overscan_subtract(utds,**kwargs)
	make_2d_biases(utds,**kwargs)
	biasMap = get_bias_map(utds)
	make_dome_flats(fileMap,biasMap,utds,**kwargs)
	if False:
		utd,filt,flatNum = '201404025','g',1
		flatFn = caldir+'flat_%s_%s_%d.fits' % (utd,filt,flatNum)
		bpMaskFile = os.path.join(caldir,'badpix_master.fits')
		build_mask_from_flat(flatFn,bpMaskFile,**kwargs)
	flatMap = get_flat_map(utds)
	# XXX propagate bpm
	process_all(fileMap,biasMap,flatMap,utds,**kwargs)
	make_supersky_flats(fileMap,utds,**kwargs)

if __name__=='__main__':
	rmpipe()

