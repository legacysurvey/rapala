#!/usr/bin/env python

import os
import re
import glob
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio

import bokutil
from bokoscan import BokOverscanSubtract
import bokproc
import badpixels

# XXX
from astrotools.idmstuff import loadpath
loadpath()
import boklog
logs = boklog.load_Bok_logs()

reduxVersion = 'bokrm_v0.1'

class RMFileNameMap(bokutil.FileNameMap):
	def __init__(self,rawDir,procDir,newSuffix=None,fromRaw=False):
		self.newSuffix = '' if newSuffix is None else newSuffix
		self.fromRaw = fromRaw
		self.rawDir = rawDir
		self.procDir = procDir
	def __call__(self,fileName):
		fn = fileName+self.newSuffix+'.fits'
		if self.fromRaw:
			fpath = os.path.join(self.rawDir,fn)
			if not os.path.exists(fpath):
				# maybe it is gzipped
				fpath = os.path.join(self.rawDir,fn+'.gz')
			return fpath
		else:
			return os.path.join(self.procDir,fn)

class FileMgr(object):
	def __init__(self,rawDir,procDir):
		self.rawDir = rawDir
		self.procDir = procDir
		self.calDir = os.path.join(self.procDir,'cals')
		self.diagDir = os.path.join(self.procDir,'diagnostics')
		self.masterBpMaskFits = None
		self.masterBpMaskFn = os.path.join(self.calDir,
		                                   'badpix_master.fits')
		self.masterBpMask4Fits = None
		self.masterBpMask4Fn = os.path.join(self.calDir,
		                                    'badpix_master_4ccd.fits')
		self.utDates = sorted(logs.keys())
		self.filt = 'gi'
	def getRawDir(self):
		return self.rawDir
	def getProcDir(self):
		return self.procDir
	def getCalDir(self):
		return self.calDir
	def getDiagDir(self):
		return self.diagDir
	def setUtDates(self,utDates):
		self.utDates = utDates
	def getUtDates(self):
		return self.utDates
	def setFilters(self,filt):
		self.filt = filt
	def getFilters(self):
		return self.filt
	def __call__(self,t,output=True):
		if t == 'raw':
			return RMFileNameMap(self.rawDir,self.procDir,fromRaw=True)
		elif t == 'MasterBadPixMask':
			if self.masterBpMaskFits is None:
				self.masterBpMaskFits = fitsio.FITS(self.masterBpMaskFn)
			return self.masterBpMaskFits
		elif t == 'MasterBadPixMask4':
			if self.masterBpMask4Fits is None:
				self.masterBpMask4Fits = fitsio.FITS(self.masterBpMask4Fn)
			return self.masterBpMask4Fits
		else:
			raise ValueError

class ProcessInPlace(FileMgr):
	def __init__(self,rawDir,procDir):
		super(ProcessInPlace,self).__init__(rawDir,procDir)
		self.fmap = RMFileNameMap(self.rawDir,self.procDir)
		# hacky to add oscan here, because files are pulled from the raw
		# directory for overscan subtraction, and need to be mapped to the
		# output directory
		self.fremap = {'oscan':'','pass1cat':'.cat1','skymask':'.skymsk',
		               '_sky':'_tmpsky'}
	def __call__(self,t,output=True):
		try:
			return super(ProcessInPlace,self).__call__(t,output)
		except ValueError:
			if t in self.fremap:
				# some files need to be remapped even when processing in-place
				return RMFileNameMap(self.rawDir,self.procDir,self.fremap[t])
			if output:
				return None
			else:
				return self.fmap

class ProcessToNewFiles(FileMgr):
	def __init__(self,rawDir,procDir):
		super(ProcessToNewFiles,self).__init__(rawDir,procDir)
		self.fmap = {'oscan':'','bias':'_b','proc':'_p','comb':'_c',
		             'pass1cat':'.cat1','skymask':'.skymsk',
		             '_sky':'_tmpsky','sky':'_s','proc2':'_q'}
	def __call__(self,t,output=True):
		try:
			return super(ProcessToNewFiles,self).__call__(t,output)
		except ValueError:
			if t == 'raw':
				return RMFileNameMap(self.rawDir,self.procDir,fromRaw=True)
			else:
				return RMFileNameMap(self.rawDir,self.procDir,self.fmap[t])

def get_utds(utds=None):
	if utds is None:
		_utds = sorted(logs.keys())
	elif type(utds) is str:
		_utds = [utds,]
	else:
		_utds = utds
	return _utds

def get_files(logs,utds=None,imType=None,filt=None,im_range=None,
              exclude_objs=None,addBiases=True):
	utds = get_utds(utds)
	files = []
	for utd in utds:
		if im_range is None:
			im_range = (0,len(logs[utd]))
		frameNum = np.arange(len(logs[utd]))
		is_range = (im_range[0] <= frameNum) & (frameNum <= im_range[1])
		if imType is None:
			is_type = True
		else:
			is_type = logs[utd]['imType'] == imType
		if filt is None:
			is_filt = True
		else:
			is_filt = logs[utd]['filter'] == filt
			if addBiases:
				# special case to include bias frames regardless of what
				# filter was in place when they were taken
				is_filt |= logs[utd]['imType'] == 'zero'
		exclude = np.zeros_like(is_range)
		if exclude_objs is not None:
			for objnm in exclude_objs:
				exclude[logs[utd]['objectName']==objnm] = True
		ii = np.where(is_range & is_type & is_filt & ~exclude)[0]
		if len(ii) > 0:
			files.append(char_add('ut'+utd+'/',logs[utd]['fileName'][ii]))
	if len(files)==0:
		return None
	else:
		return np.concatenate(files)

def get_bias_map(file_map):
	biasMap = {}
	biasPattern = os.path.join(file_map.getCalDir(),'bias_*.fits')
	biasFiles = sorted(glob.glob(biasPattern))
	bias2utd = ".*bias_(\d+).*"
	biasUtds = np.array([int(re.match(bias2utd,fn).groups()[0])
	                       for fn in biasFiles])
	for i,utd in enumerate(file_map.getUtDates()):
		j = np.argmin(np.abs(int(utd)-biasUtds))
		biasFile = biasFiles[j]
		for f in get_files(logs,utd,filt=file_map.getFilters()):
			biasMap[f] = biasFile
	return biasMap

def get_flat_map(file_map):
	flatMap = {}
	flatPattern = os.path.join(file_map.getCalDir(),'flat_????????_?_?.fits')
	flatFiles = sorted(glob.glob(flatPattern))
	flat2utdfilt = ".*flat_(\d+)_(\w+)_.*"
	utdfilt = [ re.match(flat2utdfilt,fn).groups() for fn in flatFiles]
	flatUtds = np.array([int(utd) for utd,filt in utdfilt])
	flatFilt = np.array([filt for utd,filt in utdfilt])
	for filt in file_map.getFilters():
		for i,utd in enumerate(file_map.getUtDates()):
			files = get_files(logs,utd,filt=filt)
			if files is None:
				continue
			jj = np.where(flatFilt==filt)[0]
			j = np.argmin(np.abs(int(utd)-flatUtds[jj]))
			flatFile = flatFiles[jj[j]]
			for f in files:
				flatMap[f] = flatFile
	return flatMap

def makeccd4image(file_map,inputFile,**kwargs):
	ccd4map = bokutil.FileNameMap(file_map.getCalDir(),'_4ccd')
	bokproc.combine_ccds([inputFile,],output_map=ccd4map,**kwargs)

def overscan_subtract(file_map,addBiases=True,**kwargs):
	oscanSubtract = BokOverscanSubtract(input_map=file_map('raw',False),
                                        output_map=file_map('oscan'),
                                        **kwargs)
	files = get_files(logs,file_map.getUtDates(),filt=file_map.getFilters(),
	                  addBiases=addBiases)
	oscanSubtract.process_files(files)

def make_2d_biases(file_map,nSkip=2,reject='sigma_clip',
                   writeccdim=False,**kwargs):
	biasStack = bokproc.BokBiasStack(input_map=file_map('oscan',False),
	                                 reject=reject,
                                     **kwargs)
	for utd in file_map.getUtDates():
		files = get_files(logs,utd,imType='zero')
		if files is None:
			continue
		files = files[nSkip:]
		biasNum = 1
		biasFile = os.path.join(file_map.getCalDir(),
		                        'bias_%s_%d.fits' % (utd,biasNum))
		biasStack.stack(files,biasFile)
		if writeccdim:
			makeccd4image(file_map,biasFile,**kwargs)

def make_dome_flats(file_map,bias_map,
                    nSkip=1,reject='sigma_clip',writeccdim=False,**kwargs):
	bias2Dsub = bokproc.BokCCDProcess(bias_map,
	                                  input_map=file_map('oscan',False),
	                                  output_map=file_map('bias'),
	                                  header_key='BIAS2D',
	                                  **kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=file_map('bias',False),
	                                     **kwargs)
	for utd in file_map.getUtDates():
		for filt in file_map.getFilters():
			files = get_files(logs,utd,imType='flat',filt=filt)
			if files is None:
				continue
			files = files[nSkip:]
			bias2Dsub.process_files(files)
			flatNum = 1
			flatFile = os.path.join(file_map.getCalDir(),
			                        'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
			flatStack.stack(files,flatFile)
			if writeccdim:
				makeccd4image(file_map,flatFile,**kwargs)

def make_bad_pixel_masks(file_map,**kwargs):
	utd,filt,flatNum = '20140425','g',1
	flatFn = os.path.join(file_map.getCalDir(),
	                      'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
	bpMaskFile = os.path.join(file_map.getCalDir(),'badpix_master.fits')
	badpixels.build_mask_from_flat(flatFn,bpMaskFile,
	             normed_flat_file=flatFn.replace('.fits','_normed.fits'),
	             normed_flat_fit_file=flatFn.replace('.fits','_fit.fits'),
	             binned_flat_file=flatFn.replace('.fits','_binned.fits'),
	             )#,**kwargs)
	makeccd4image(file_map,bpMaskFile,**kwargs)

def balance_gains(file_map,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=file_map('proc',False),
	                                     mask_map=file_map('MasterBadPixMask'),
	                                                **kwargs)
	gainMap = {'corrections':{},'skyvals':{}}
	for utd in file_map.getUtDates():
		for filt in file_map.getFilters():
			files = get_files(logs,utd,imType='object',filt=filt)
			if files is None:
				continue
			diagfile = os.path.join(fileMap.getDiagDir(),
			                        'gainbal_%s_%s.npz'%(utd,filt))
			if os.path.exists(diagfile):
				gainDat = np.load(diagfile)
				gainCor = gainDat['gainCor']
				skyV = gainDat['skys']
			else:
				gainBalance.process_files(files)
				gainCor = gainBalance.calc_mean_corrections()
				gainCorV,skyV = gainBalance.get_values()
			for f,skyv in zip(files,skyV):
				gainMap['corrections'][f] = gainCor
				gainMap['skyvals'][f] = skyv
			gainBalance.reset()
			if not os.path.exists(diagfile):
				np.savez(diagfile,gains=gainCorV,skys=skyV,gainCor=gainCor)
	return gainMap

def process_all(file_map,bias_map,flat_map,
                fixpix=False,nocombine=False,**kwargs):
	# 1. basic processing (bias and flat-field correction, fixpix, 
	#    nominal gain correction
	proc = bokproc.BokCCDProcess(bias_map,
	                             flat_map,
	                             input_map=file_map('oscan',False),
	                             output_map=file_map('proc'),
	                             mask_map=file_map('MasterBadPixMask'),
	                             fixpix=fixpix,
	                             **kwargs)
	files = get_files(logs,file_map.getUtDates(),filt=file_map.getFilters(),
	                  imType='object')
	proc.process_files(files)
	if nocombine:
		return
	# 2. balance gains using background counts
	gainMap = balance_gains(file_map,**kwargs)
	# 3. combine per-amp images (16) into CCD images (4)
	bokproc.combine_ccds(files,
	                     input_map=file_map('proc',False),
	                     output_map=file_map('comb'),
	                     gain_map=gainMap,
	                     **kwargs)

def make_supersky_flats(file_map,skysub=True,**kwargs):
	if skysub:
		skySub = bokproc.BokSkySubtract(input_map=file_map('comb',False),
		                                output_map=file_map('_sky'),
		                                mask_map=file_map('skymask'))
		skySub.add_mask(file_map('MasterBadPixMask4'))
		stackin = file_map('_sky',False)
	else:
		stackin = file_map('comb',False)
	skyFlatMask = bokproc.BokGenerateSkyFlatMasks(
	                                    input_map=file_map('comb',False),
	                                    output_map=file_map('skymask'),
	                                    mask_map=file_map('MasterBadPixMask4'))
	caldir = file_map.getCalDir()
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=file_map('skymask'),
	                    exposure_time_map=bokutil.FileNameMap(caldir,'.exp'),
	                       raw_stack_file=bokutil.FileNameMap(caldir,'_raw'),
	                                        header_bad_key='BADSKY')
	skyFlatStack.set_badpixelmask(file_map('MasterBadPixMask4'))
	for utd in file_map.getUtDates():
		for filt in file_map.getFilters():
			# exclude RM10 and RM11 because they are swamped by a bright star
			files = get_files(logs,utd,imType='object',filt=filt,
			                  )#exclude_objs=['rm10','rm11'])
			if files is None:
				continue
			# XXX need to use the bad pix mask for weighting
			#bokproc.sextract_pass1(files,
			#                       input_map=file_map('comb',False),
			#                       catalog_map=file_map('pass1cat'),
			#                       object_mask_map=file_map('skymask'),
			#                       **kwargs)
			skyFlatMask.process_files(files)
			if skysub:
				skySub.process_files(files)
			outfn = os.path.join(file_map.getCalDir(),
			                     'skyflat_%s_%s.fits' % (utd,filt))
			skyFlatStack.stack(files,outfn)

# process round 2: illum corr and sky sub

def make_images():
	import matplotlib.pyplot as plt
	files = get_files(logs,'20140427',imType='object',filt='g')
	_fmap = RMFileNameMap()
	plt.ioff()
	for ff in files:
		f = _fmap(ff)
		if not os.path.exists(f.replace('.fits','_tmpsky.fits')):
			continue
		if os.path.exists(f.replace('.fits','.png')):
			continue
		bokproc.make_fov_image_fromfile(f.replace('.fits','_tmpsky.fits'),
		                                f.replace('.fits','.png'),
		                                mask=f.replace('.fits','.skymsk.fits'))
	plt.ion()

all_process_steps = ['oscan','bias2d','flat2d','bpmask',
                     'proc1','skyflat','proc2']

def create_file_map(rawDir,procDir,utds,bands,newfiles):
	# set default file paths
	if rawDir is None:
		rawDir = os.environ['BOK90PRIMERAWDIR']
	if procDir is None:
		procDir = os.path.join(os.environ['BOK90PRIMEOUTDIR'],reduxVersion)
	# create the file manager object
	if newfiles:
		fileMap = ProcessToNewFiles(rawDir,procDir) 
	else:
		fileMap = ProcessInPlace(rawDir,procDir) 
	# restrict processing to specified UT dates and filters
	if utds is not None:
		fileMap.setUtDates(utds)
	if bands is not None:
		fileMap.setFilters(bands)
	# create output directories for processed files
	for d in [procDir,fileMap.getCalDir(),fileMap.getDiagDir()]:
		if not os.path.exists(d):
			os.mkdir(d)
	for utd in utds:
		utdir = os.path.join(procDir,'ut'+utd)
		if not os.path.exists(utdir): os.mkdir(utdir)
	return fileMap

def rmpipe(fileMap,redo,steps,verbose,**kwargs):
	pipekwargs = {'clobber':redo,'verbose':verbose}
	# fixpix is sticking nan's into the images in unmasked pixels (???)
	fixpix = False #True
	timerLog = bokutil.TimerLog()
	if 'oscan' in steps:
		overscan_subtract(fileMap,addBiases=True,**pipekwargs)
		timerLog('overscans')
	if 'bias2d' in steps:
		make_2d_biases(fileMap,writeccdim=True,**pipekwargs)
		timerLog('2d biases')
	biasMap = get_bias_map(fileMap)
	if 'flat2d' in steps:
		make_dome_flats(fileMap,biasMap,writeccdim=True,**pipekwargs)
		timerLog('dome flats')
	if 'bpmask' in steps:
		make_bad_pixel_masks(fileMap)
		timerLog('bad pixel masks')
	if kwargs.get('noflatcorr',False):
		flatMap = None
	else:
		flatMap = get_flat_map(fileMap)
	if 'proc1' in steps:
		process_all(fileMap,biasMap,flatMap,
		            fixpix=fixpix,nocombine=kwargs.get('nocombine'),
		            **pipekwargs)
		timerLog('ccdproc')
	if 'skyflat' in steps:
		make_supersky_flats(fileMap,**pipekwargs)
		timerLog('supersky flats')
	if 'proc2' in steps:
		# XXX for testing
		#fileMap = ProcessToNewFiles()
		timerLog('process2')
	timerLog.dump()

if __name__=='__main__':
	import sys
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--band',type=str,default=None,
	                help='band to process (g or i) [default=both]')
	parser.add_argument('-i','--images',action='store_true',
	                help='make png images [default=False]')
	parser.add_argument('-n','--newfiles',action='store_true',
	                help='process to new files (not in-place)')
	parser.add_argument('-o','--output',type=str,default=None,
	                help='output directory [default=$BOK90PRIMEOUTDIR]')
	parser.add_argument('-r','--rawdir',type=str,default=None,
	                help='raw data directory [default=$BOK90PRIMERAWDIR]')
	parser.add_argument('-R','--redo',action='store_true',
	                help='redo (overwrite existing files)')
	parser.add_argument('-s','--steps',type=str,default=None,
	                help='processing steps to execute [default=all]')
	parser.add_argument('-S','--stepto',type=str,default=None,
	                help='process until this step [default=last]')
	parser.add_argument('-u','--utdate',type=str,default=None,
	                help='UT date(s) to process [default=all]')
	parser.add_argument('-v','--verbose',action='count',
	                help='increase output verbosity')
	parser.add_argument('--noflatcorr',action='store_true',
	                help='do not apply flat correction')
	parser.add_argument('--nocombine',action='store_true',
	                help='do not combine into CCD images')
	args = parser.parse_args()
	if args.utdate is None:
		utds = None
	else:
		utds = args.utdate.split(',')
	if args.steps is None:
		if args.stepto is None:
			steps = all_process_steps
		else:
			steps = all_process_steps[:all_process_steps.index(args.stepto)+1]
	else:
		steps = args.steps.split(',')
	verbose = 0 if args.verbose is None else args.verbose
	fileMap = create_file_map(args.rawdir,args.output,
	                          utds,args.band,args.newfiles)
	if args.images:
		makeimages()
	else:
		rmpipe(fileMap,args.redo,steps,verbose,
		       noflatcorr=args.noflatcorr,
		       nocombine=args.nocombine)

