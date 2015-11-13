#!/usr/bin/env python

import os
import re
import glob
from copy import copy
import multiprocessing
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio

from bokpipe import *
from bokpipe import __version__ as pipeVersion

# XXX
from astrotools.idmstuff import loadpath
loadpath()
import boklog
logs = boklog.load_Bok_logs()

all_process_steps = ['oscan','bias2d','flat2d','bpmask',
                     'proc1','skyflat','proc2']

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
		self.masterBpMaskFn = 'badpix_master.fits'
		self.masterBpMaskFits = None
		self.masterBpMask4Fn = 'badpix_master_4ccd.fits'
		self.masterBpMask4Fits = None
		self.masterRampCorrFn = None
		self.masterRampCorrFits = None
		self.utDates = sorted(logs.keys())
		self.filt = 'gi'
		self.frames = None
		self._curUtDate = None
		self._curFilt = None
	def getRawDir(self):
		return self.rawDir
	def getProcDir(self):
		return self.procDir
	def getCalDir(self):
		return self.calDir
	def setCalDir(self,calDir):
		self.calDir = calDir
	def getDiagDir(self):
		return self.diagDir
	def setDiagDir(self,diagDir):
		self.diagDir = diagDir
	def setUtDates(self,utDates):
		self.utDates = utDates
	def getUtDates(self):
		return self.utDates
	def iterUtDates(self):
		for u in self.utDates:
			self._curUtDate = u
			yield u
		self._curUtDate = None
	def setFilters(self,filt):
		self.filt = filt
	def getFilters(self):
		return self.filt
	def iterFilters(self):
		for f in self.filt:
			self._curFilt = f
			yield f
		self._curFilt = None
	def setFrames(self,frames):
		self.frames = frames
	def setRampCorrFile(self,rampCorrFn):
		self.masterRampCorrFn = rampCorrFn
	def __call__(self,t,output=True):
		if t == 'raw':
			return RMFileNameMap(self.rawDir,self.procDir,fromRaw=True)
		elif t == 'MasterBadPixMask':
			if self.masterBpMaskFits is None:
				fn = os.path.join(self.calDir,self.masterBpMaskFn)
				self.masterBpMaskFits = fitsio.FITS(fn)
			return self.masterBpMaskFits
		elif t == 'MasterBadPixMask4':
			if self.masterBpMask4Fits is None:
				fn = os.path.join(self.calDir,self.masterBpMask4Fn)
				self.masterBpMask4Fits = fitsio.FITS(fn)
			return self.masterBpMask4Fits
		elif t == 'BiasRampCorrection':
			if self.masterRampCorrFits is None \
			      and self.masterRampCorrFn is not None:
				self.masterRampCorrFits = fitsio.FITS(self.masterRampCorrFn)
			return self.masterRampCorrFits
		else:
			raise ValueError
	def getFiles(self,imType=None,utd=None,filt=None,
	             im_range=None,exclude_objs=None,with_frames=False):
		if utd is not None:
			utds = [utd] if type(utd) is str else utd
		elif self._curUtDate is None:
			utds = self.getUtDates() 
		else:
			utds = [self._curUtDate]
		files,frames = [],[]
		for utd in utds:
			nFrames = len(logs[utd])
			if im_range is None and self.frames is None:
				# need one boolean array to be full length
				is_range = np.ones(nFrames,dtype=bool)
			else:
				if im_range is not None:
					i1,i2 = im_range
				else:
					i1,i2 = self.frames
				frameNum = np.arange(nFrames)
				is_range = ( (i1 <= frameNum) & (frameNum <= i2) )
			if imType is None:
				is_type = True
			else:
				is_type = logs[utd]['imType'] == imType
			if filt is None and self._curFilt is None and self.filt=='gi':
				# using all filters (good thing there aren't >2!)
				is_filt = True
			else:
				# restricted to one filter
				if filt is not None:
					f = filt
				elif self._curFilt is not None:
					f = self._curFilt
				else:
					f = self.filt
				is_filt = logs[utd]['filter'] == f
				if imType is None:
					# special case to include bias frames regardless of 
					# what filter was in place when they were taken
					is_filt |= logs[utd]['imType'] == 'zero'
			if exclude_objs is None:
				exclude = False
			else:
				exclude = np.zeros(nFrames,dtype=bool)
				for objnm in exclude_objs:
					exclude[logs[utd]['objectName']==objnm] = True
			ii = np.where(is_range & is_type & is_filt & ~exclude)[0]
			if len(ii) > 0:
				files.append(char_add('ut'+utd+'/',logs[utd]['fileName'][ii]))
				if with_frames:
					frames.append(ii)
				# XXX could use utds here for getting nearest cal
		if len(files)==0:
			files = None
			if with_frames:
				frames = None
		else:
			files = np.concatenate(files)
			if with_frames:
				frames = np.concatenate(frames)
		if with_frames:
			return files,frames
		else:
			return files

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

def get_bias_map(file_map):
	biasMap = {}
	biasPattern = os.path.join(file_map.getCalDir(),'bias_*.fits')
	biasFiles = sorted(glob.glob(biasPattern))
	bias2utd = ".*bias_(\d+).*"
	biasUtds = np.array([int(re.match(bias2utd,fn).groups()[0])
	                       for fn in biasFiles])
	for utd in file_map.iterUtDates():
		j = np.argmin(np.abs(int(utd)-biasUtds))
		biasFile = biasFiles[j]
		for f in file_map.getFiles():
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
	for filt in file_map.iterFilters():
		for utd in file_map.iterUtDates():
			files = file_map.getFiles()
			if files is None:
				continue
			jj = np.where(flatFilt==filt)[0]
			try:
				j = np.argmin(np.abs(int(utd)-flatUtds[jj]))
			except ValueError:
				raise ValueError('no match for %s:%s in %s:s' % 
				                 (utd,filt,flatUtds,flatFilt))
			flatFile = flatFiles[jj[j]]
			for f in files:
				flatMap[f] = flatFile
	return flatMap

def makeccd4image(file_map,inputFile,**kwargs):
	ccd4map = bokutil.FileNameMap(file_map.getCalDir(),'_4ccd')
	bokproc.combine_ccds([inputFile,],output_map=ccd4map,**kwargs)

def overscan_subtract(file_map,**kwargs):
	oscanSubtract = BokOverscanSubtract(input_map=file_map('raw',False),
                                        output_map=file_map('oscan'),
                                        **kwargs)
	oscanSubtract.process_files(file_map.getFiles())

def make_2d_biases(file_map,nSkip=2,reject='sigma_clip',
                   writeccdim=False,**kwargs):
	biasStack = bokproc.BokBiasStack(input_map=file_map('oscan',False),
	                                 reject=reject,
                                     **kwargs)
	for utd in file_map.iterUtDates():
		files,frames = file_map.getFiles(imType='zero',with_frames=True)
		if files is None:
			continue
		splits = np.where(np.diff(frames)>1)[0]
		if len(splits)==0:
			bias_seqs = [ files ]
		else:
			bias_seqs = np.split(files,splits+1)
		for biasNum,biasFiles in enumerate(bias_seqs,start=1):
			if len(biasFiles) < 5: # XXX hardcoded
				continue
			biasFile = os.path.join(file_map.getCalDir(),
			                        'bias_%s_%d.fits' % (utd,biasNum))
			biasStack.stack(biasFiles[nSkip:],biasFile)
			if writeccdim:
				makeccd4image(file_map,biasFile,**kwargs)

def make_dome_flats(file_map,bias_map,
                    nSkip=1,reject='sigma_clip',writeccdim=False,
	                usepixflat=True,**kwargs):
	bias2Dsub = bokproc.BokCCDProcess(bias_map,
	                                  input_map=file_map('oscan',False),
	                                  output_map=file_map('bias'),
	                                  header_key='BIAS2D',
	                                  **kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=file_map('bias',False),
	                                     **kwargs)
	if usepixflat:
		normFlat = bokproc.NormalizeFlat(**kwargs)
	for utd in file_map.iterUtDates():
		for filt in file_map.iterFilters():
			files,frames = file_map.getFiles(imType='flat',with_frames=True)
			if files is None:
				continue
			splits = np.where(np.diff(frames)>1)[0]
			if len(splits)==0:
				flat_seqs = [ files ]
			else:
				flat_seqs = np.split(files,splits+1)
			for flatNum,flatFiles in enumerate(flat_seqs,start=1):
				flatFile = os.path.join(file_map.getCalDir(),
			                        'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
				bias2Dsub.process_files(flatFiles)
				flatStack.stack(flatFiles[nSkip:],flatFile)
				if usepixflat:
					normFlat.process_files([flatFile])
				if writeccdim:
					makeccd4image(file_map,flatFile,**kwargs)

def make_bad_pixel_masks(file_map,**kwargs):
	# XXX hardcoded
	utd,filt,flatNum = '20140425','g',1
	flatFn = os.path.join(file_map.getCalDir(),
	                      'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
	bpMaskFile = os.path.join(file_map.getCalDir(),'badpix_master.fits')
	badpixels.build_mask_from_flat(flatFn,bpMaskFile)#,**kwargs)
	makeccd4image(file_map,bpMaskFile,**kwargs)

def balance_gains(file_map,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=file_map('proc',False),
	                                     mask_map=file_map('MasterBadPixMask'),
	                                                **kwargs)
	gainMap = {'corrections':{},'skyvals':{}}
	for utd in file_map.iterUtDates():
		for filt in file_map.iterFilters():
			files = file_map.getFiles(imType='object')
			if files is None:
				continue
			diagfile = os.path.join(file_map.getDiagDir(),
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
	                             ramp_map=file_map('BiasRampCorrection'),
	                             fixpix=fixpix,
	                             **kwargs)
	files = file_map.getFiles(imType='object')
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
	for utd in file_map.iterUtDates():
		for filt in file_map.iterFilters():
			# exclude RM10 and RM11 because they are swamped by a bright star
			files = file_map.getFiles(imType='object',
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

def make_images(file_map,imtype='comb',msktype=None):
	import matplotlib.pyplot as plt
	files = file_map.getFiles(imType='object')
	_fmap = file_map(imtype)
	if msktype=='badpix':
		msktype = 'MasterBadPixMask4'
	if msktype==None:
		maskmap = lambda f: None
	else:
		maskmap = file_map(msktype)
	imdir = os.path.join(file_map.getProcDir(),'images')
	if not os.path.exists(imdir):
		os.mkdir(imdir)
	plt.ioff()
	for ff in files:
		f = _fmap(ff)
		if not os.path.exists(f):
			continue
		imgfile = os.path.basename(f).replace('.fits','.png')
		imgfile = os.path.join(imdir,imgfile)
		if os.path.exists(imgfile):
			continue
		bokproc.make_fov_image_fromfile(f,imgfile,mask=maskmap(ff))
	plt.ion()

def create_file_map(rawDir,procDir,utds,bands,newfiles):
	# set default file paths
	if rawDir is None:
		rawDir = os.environ['BOK90PRIMERAWDIR']
	if procDir is None:
		procDir = os.path.join(os.environ['BOK90PRIMEOUTDIR'],pipeVersion)
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
	for utd in fileMap.getUtDates():
		utdir = os.path.join(procDir,'ut'+utd)
		if not os.path.exists(utdir): os.mkdir(utdir)
	return fileMap

def rmpipe(fileMap,redo,steps,verbose,**kwargs):
	pipekwargs = {'clobber':redo,'verbose':verbose}
	# fixpix is sticking nan's into the images in unmasked pixels (???)
	fixpix = False #True
	writeccdims = kwargs.get('calccdims',False)
	timerLog = bokutil.TimerLog()
	biasMap = None
	if 'oscan' in steps:
		overscan_subtract(fileMap,**pipekwargs)
		timerLog('overscans')
	if 'bias2d' in steps:
		make_2d_biases(fileMap,writeccdim=writeccdims,**pipekwargs)
		timerLog('2d biases')
	if 'flat2d' in steps:
		biasMap = get_bias_map(fileMap)
		make_dome_flats(fileMap,biasMap,writeccdim=writeccdims,
		                usepixflat=not kwargs.get('nousepixflat',False),
		                **pipekwargs)
		timerLog('dome flats')
	if 'bpmask' in steps:
		make_bad_pixel_masks(fileMap)
		timerLog('bad pixel masks')
	if 'proc1' in steps:
		if biasMap is None:
			biasMap = get_bias_map(fileMap)
		if kwargs.get('noflatcorr',False):
			flatMap = None
		else:
			flatMap = get_flat_map(fileMap)
		if not kwargs.get('norampcorr',False):
			fileMap.setRampCorrFile(os.path.join(fileMap.getCalDir(),
			                                     'biasramp.fits'))
		process_all(fileMap,biasMap,flatMap,
		            fixpix=fixpix,
		            nocombine=kwargs.get('nocombine'),
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

def rmpipe_poormp(nProc,fileMap,*args,**kwargs):
	def chunks(l, n):
		nstep = int(round(len(l)/float(n)))
		for i in xrange(0, len(l), nstep):
			yield l[i:i+nstep]
	utdSets = chunks(fileMap.getUtDates(),nProc)
	jobs = []
	for i,utds in enumerate(utdSets):
		fmap = copy(fileMap)
		fmap.setUtDates(utds)
		_args = (fmap,) + args
		p = multiprocessing.Process(target=rmpipe,args=_args,kwargs=kwargs)
		jobs.append(p)
		p.start()

if __name__=='__main__':
	import sys
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--band',type=str,default=None,
	                help='band to process (g or i) [default=both]')
	parser.add_argument('-c','--caldir',type=str,default=None,
	                help='set calibration directory')
	parser.add_argument('-f','--frames',type=str,default=None,
	                help='frames to process (i1,i2) [default=all]')
	parser.add_argument('-i','--images',type=str,default=None,
	                help='make png images (imtype,[msktype]) [default=no]')
	parser.add_argument('-n','--newfiles',action='store_true',
	                help='process to new files (not in-place)')
	parser.add_argument('-o','--output',type=str,default=None,
	                help='output directory [default=$BOK90PRIMEOUTDIR]')
	parser.add_argument('-p','--processes',type=int,default=1,
	                help='number of processes to use [default=single]')
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
	parser.add_argument('--calccdims',action='store_true',
	                help='generate CCD-combined images for calibration data')
	parser.add_argument('--noflatcorr',action='store_true',
	                help='do not apply flat correction')
	parser.add_argument('--nocombine',action='store_true',
	                help='do not combine into CCD images')
	parser.add_argument('--norampcorr',action='store_true',
	                help='do not attempt to correct bias ramp')
	parser.add_argument('--nousepixflat',action='store_true',
	                help='do not use normalized pixel flat')
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
	if args.frames is not None:
		fileMap.setFrames(tuple([int(_f) for _f in args.frames.split(',')]))
	if args.caldir is not None:
		fileMap.setCalDir(os.path.join(args.caldir,'cals'))
		fileMap.setDiagDir(os.path.join(args.caldir,'diagnostics'))
	if args.images is not None:
		make_images(fileMap,*args.images.split(','))
	elif args.processes > 1:
		rmpipe_poormp(args.processes,
		              fileMap,args.redo,steps,verbose,
		              noflatcorr=args.noflatcorr,
		              nocombine=args.nocombine,
		              calccdims=args.calccdims,
		              norampcorr=args.norampcorr,
		              nousepixflat=args.nousepixflat)
	else:
		rmpipe(fileMap,args.redo,steps,verbose,
		       noflatcorr=args.noflatcorr,
		       nocombine=args.nocombine,
		       calccdims=args.calccdims,
		       norampcorr=args.norampcorr,
		       nousepixflat=args.nousepixflat)

