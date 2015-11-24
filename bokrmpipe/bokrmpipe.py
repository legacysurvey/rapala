#!/usr/bin/env python

import os
import re
import glob
from copy import copy
import multiprocessing
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio
from astropy.table import Table

from bokpipe import *
from bokpipe import __version__ as pipeVersion

import bokrampcorr
import bokillumcorr

all_process_steps = ['oscan','bias2d','flat2d','bpmask',
                     'proc1','skyflat','proc2','wcs','cat']

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
	def __init__(self,obsDb,rawDir,procDir):
		self.rawDir = rawDir
		self.procDir = procDir
		self.calDir = os.path.join(self.procDir,'cals')
		self.diagDir = os.path.join(self.procDir,'diagnostics')
		self.masterBpMaskFn = 'badpix_master.fits'
		self.masterBpMaskFits = None
		self.masterBpMask4Fn = 'badpix_master_4ccd.fits'
		self.masterBpMask4Fits = None
		self.masterRampCorrFn = 'biasramp.fits'
		self.masterRampCorrFits = None
		self.illumCorrFn = 'illumination.fits'
		self.illumCorrFits = None
		self.skyFlatFn = 'skyflat_g.fits'
		self.skyFlatFits = None
		self.obsDb = obsDb
		self.allUtDates = np.unique(self.obsDb['utDate'])
		self.utDates = self.allUtDates
		self.filt = 'gi'
		self.frames = None
		self.frameList = None
		self.imType = None
		self._curUtDate = None
		self._curFilt = None
		self._tmpInput = False
		self._tmpOutput = False
		self._tmpDir = os.path.join(self.procDir,'tmp')
	def getRawDir(self):
		return self.rawDir
	def getProcDir(self):
		return self.procDir
	def setTmpInput(self):
		self._tmpInput = True
	def setTmpOutput(self):
		self._tmpOutput = True
		if not os.path.exists(self._tmpDir):
			os.mkdir(self._tmpDir)
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
	def setFile(self,fileName,utDate=None):
		isUtd = True if utDate is None else self.obsDb['utDate']==utDate
		self.frameList = np.where( isUtd & 
		                           (self.obsDb['fileName']==fileName) )[0]
	def setFrameList(self,utDates,fileNames):
		frames = [ np.where( (self.obsDb['utDate']==utd) &
		                     (self.obsDb['fileName']==f) )[0][0]
		             for utd,f in zip(utDates,fileNames) ]
		self.frameList = np.array(frames)
	def setImageType(self,imType):
		self.imType = imType
	def setScampRefCatDir(self,refCatDir):
		self.refCatDir = refCatDir
		if not os.path.exists(self.refCatDir):
			os.mkdir(self.refCatDir)
	def getScampRefCat(self,fieldName):
		refCatFn = 'scampref_%s.cat' % fieldName
		return os.path.join(self.refCatDir,refCatFn)
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
			if self.masterRampCorrFits is None:
				fn = os.path.join(self.calDir,self.masterRampCorrFn)
				self.masterRampCorrFits = fitsio.FITS(fn)
			return self.masterRampCorrFits
		elif t == 'IllumCorrImage':
			if self.illumCorrFits is None:
				fn = os.path.join(self.calDir,self.illumCorrFn)
				self.illumCorrFits = fitsio.FITS(fn)
			return self.illumCorrFits
		elif t == 'DarkSkyFlatImage':
			if self.skyFlatFits is None:
				fn = os.path.join(self.calDir,self.skyFlatFn)
				self.skyFlatFits = fitsio.FITS(fn)
			return self.skyFlatFits
		else:
			raise ValueError
	def getFiles(self,imType=None,utd=None,filt=None,
	             im_range=None,exclude_objs=None,
	             with_objnames=False,with_frames=False):
		file_sel = np.zeros(len(self.obsDb),dtype=bool)
		# select on UT date(s)
		if utd is not None:
			utds = [utd] if type(utd) is str else utd
		elif self._curUtDate is not None:
			utds = [self._curUtDate]
		else:
			utds = self.getUtDates() 
		if np.array_equal(utds,self.allUtDates):
			# all utds are valid
			file_sel[:] = True
		else:
			for utd in utds:
				file_sel |= self.obsDb['utDate'] == utd
		# restrict on image type
		if imType is not None:
			file_sel &= self.obsDb['imType'] == imType
		elif self.imType is not None:
			file_sel &= self.obsDb['imType'] == self.imType
		# restrict on filter
		if filt is not None:
			f = filt
		elif self._curFilt is not None:
			f = self._curFilt
		else:
			f = self.filt
		if f != 'gi':
			isFilt = self.obsDb['filter'] == f
			if imType is None:
				# special case to include bias frames regardless of 
				# what filter was in place when they were taken
				isFilt |= self.obsDb['imType'] == 'zero'
			file_sel &= isFilt
		# restrict on specified range of frames
		if im_range is not None:
			file_sel &= ( (self.obsDb['frameIndex'] >= im_range[0]) & 
			              (self.obsDb['frameIndex'] <= im_range[1]) )
		elif self.frames is not None:
			file_sel &= ( (self.obsDb['frameIndex'] >= self.frames[0]) & 
			              (self.obsDb['frameIndex'] <= self.frames[1]) )
		# list of objects to exclude
		if exclude_objs is not None:
			for objnm in exclude_objs:
				file_sel ^= self.obsDb['objName'] == objnm
		# finally check input frame list
		if self.frameList is not None:
			isFrame = np.zeros(len(self.obsDb),dtype=bool)
			isFrame[self.frameList] = True
			file_sel &= isFrame
		# construct the output file list
		ii = np.where(file_sel)[0]
		if len(ii) > 0:
			# XXX os.path.join for arrays?
			files = char_add(char_add(self.obsDb['utDir'][ii],'/'),
			                 self.obsDb['fileName'][ii])
			rv = [files]
			if with_objnames:
				rv.append(self.obsDb['objName'][ii])
			if with_frames:
				rv.append(ii)
			if len(rv)>1:
				return tuple(rv)
			else:
				return rv[0]
		else:
			if with_frames:
				return None,None
			else:
				return None

class ProcessInPlace(FileMgr):
	def __init__(self,obsDb,rawDir,procDir):
		super(ProcessInPlace,self).__init__(obsDb,rawDir,procDir)
		self.fmap = RMFileNameMap(self.rawDir,self.procDir)
		# hacky to add oscan here, because files are pulled from the raw
		# directory for overscan subtraction, and need to be mapped to the
		# output directory
		self.fremap = {'oscan':'','pass1cat':'.cat1',
		               'skymask':'.skymsk','skyfit':'.sky',
		               'wcscat':'.wcscat','cat':'.cat','psf':'.psf'}
	def __call__(self,t,output=True):
		if output:
			outDir = self.procDir if not self._tmpOutput else self._tmpDir
		else:
			outDir = self.procDir if not self._tmpInput else self._tmpDir
		try:
			return super(ProcessInPlace,self).__call__(t,output)
		except ValueError:
			if t in self.fremap:
				# some files need to be remapped even when processing in-place
				return RMFileNameMap(self.rawDir,outDir,self.fremap[t])
			elif output:
				return RMFileNameMap(self.rawDir,outDir)
			else:
				if self._tmpInput:
					return RMFileNameMap(self.rawDir,outDir)
				else:
					return self.fmap

class ProcessToNewFiles(FileMgr):
	def __init__(self,obsDb,rawDir,procDir):
		super(ProcessToNewFiles,self).__init__(obsDb,rawDir,procDir)
		self.fmap = {'oscan':'','bias':'_b','proc':'_p','comb':'_c',
		             'pass1cat':'.cat1','skymask':'.skymsk','skyfit':'.sky',
		             'sky':'_s','proc2':'_q','wcscat':'.wcscat',
		             'cat':'.cat','psf':'.psf'}
	def __call__(self,t,output=True):
		if output:
			outDir = self.procDir if not self._tmpOutput else self._tmpDir
		else:
			outDir = self.procDir if not self._tmpInput else self._tmpDir
		try:
			return super(ProcessToNewFiles,self).__call__(t,output)
		except ValueError:
			if t == 'raw':
				return RMFileNameMap(self.rawDir,outDir,fromRaw=True)
			else:
				return RMFileNameMap(self.rawDir,outDir,self.fmap[t])

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
		files = file_map.getFiles()
		if files is None:
			continue
		for f in files:
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
				if len(flatFiles) < 5: # XXX hardcoded
					continue
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
	build_mask_from_flat(flatFn,bpMaskFile)#,**kwargs)
	makeccd4image(file_map,bpMaskFile,**kwargs)

def balance_gains(file_map,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=file_map('proc'),
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
			if not os.path.exists(diagfile) and \
			     not kwargs.get('nosavegain',False):
				np.savez(diagfile,gains=gainCorV,skys=skyV,gainCor=gainCor)
	return gainMap

def process_all(file_map,bias_map,flat_map,
                fixpix=False,norampcorr=False,
                nocombine=False,prockey='CCDPROC',**kwargs):
	# 1. basic processing (bias and flat-field correction, fixpix, 
	#    nominal gain correction
	ramp = None if norampcorr else file_map('BiasRampCorrection')
	proc = bokproc.BokCCDProcess(input_map=file_map('oscan',False),
	                             output_map=file_map('proc'),
	                             mask_map=file_map('MasterBadPixMask'),
	                             header_key=prockey,
	                             bias=bias_map,flat=flat_map,
	                             ramp=ramp,illum=None,darksky=None,
	                             fixpix=fixpix,**kwargs)
	files = file_map.getFiles(imType='object')
	if files is None:
		return
	proc.process_files(files)
	if nocombine:
		return
	# 2. balance gains using background counts
	gainMap = balance_gains(file_map,**kwargs)
	# 3. combine per-amp images (16) into CCD images (4)
	bokproc.combine_ccds(files,
	                     input_map=file_map('proc'), # using output from above
	                     output_map=file_map('comb'),
	                     gain_map=gainMap,
	                     **kwargs)

def make_supersky_flats(file_map,**kwargs):
	caldir = file_map.getCalDir()
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=file_map('skymask'),
	                    exposure_time_map=bokutil.FileNameMap(caldir,'.exp'),
	                       raw_stack_file=bokutil.FileNameMap(caldir,'_raw'),
	                                        header_bad_key='BADSKY')
	skyFlatStack.set_badpixelmask(file_map('MasterBadPixMask4'))
	for filt in file_map.iterFilters():
		if True:
			files = file_map.getFiles(imType='object')
			if files is not None:
				outfn = os.path.join(file_map.getCalDir(),
				                     'skyflat_%s.fits' % (filt))
				skyFlatStack.stack(files,outfn)
		if False:
			# this does a sky flat each night
			for utd in file_map.iterUtDates():
				files = file_map.getFiles(imType='object')
				if files is None:
					continue
				outfn = os.path.join(file_map.getCalDir(),
				                     'skyflat_%s_%s.fits' % (utd,filt))
				skyFlatStack.stack(files,outfn)

def process_all2(file_map,skyArgs,noillumcorr=False,nodarkskycorr=False,
                 noskysub=False,prockey='CCDPRO2',save_sky=False,
                 **kwargs):
	#
	# Second round flat-field corrections
	#
	illum = None if noillumcorr else file_map('IllumCorrImage')
	darksky = None if nodarkskycorr else file_map('DarkSkyFlatImage')
	proc = bokproc.BokCCDProcess(input_map=file_map('comb',False),
	                             output_map=file_map('proc2'),
	                             mask_map=file_map('MasterBadPixMask4'),
	                             header_key=prockey,
	                             gain_multiply=False,bias=None,flat=None,
	                             ramp=None,illum=illum,darksky=darksky,
	                             fixpix=False,**kwargs)
	files = file_map.getFiles(imType='object')
	if files is None:
		return
	proc.process_files(files)
	if noskysub:
		return
	#
	# Sky subtraction
	#
	# Generate sky masks by agressively masking objects
	# old way using sextractor masks
	#bokproc.sextract_pass1(files,
	#                       input_map=file_map('comb',False),
	#                       catalog_map=file_map('pass1cat'),
	#                       object_mask_map=file_map('skymask'),
	#                       **kwargs)
	skyFlatMask = bokproc.BokGenerateSkyFlatMasks(
	                                    input_map=file_map('proc2'),
	                                    output_map=file_map('skymask'),
	                                    mask_map=file_map('MasterBadPixMask4'))
	files = file_map.getFiles(imType='object')
	skyFlatMask.process_files(files)
	skyfitmap = file_map('skyfit') if save_sky else None
	skySub = bokproc.BokSkySubtract(input_map=file_map('proc2'),
	                                output_map=file_map('sky'),
	                                mask_map=file_map('skymask'),
	                                skyfit_map=skyfitmap,**skyArgs)
	skySub.add_mask(file_map('MasterBadPixMask4'))
	stackin = file_map('sky',False)
	skySub.process_files(files)

def set_wcs(file_map,inputType='sky',keepwcscat=True,**kwargs):
	filesAndFields = file_map.getFiles(imType='object',with_objnames=True)
	for imFile,fieldName in zip(*filesAndFields):
		imageFile = file_map(inputType)(imFile)
		catFile = file_map('wcscat',output=True)(imFile)
		bokphot.sextract(imageFile,catFile,**kwargs)
		bokastrom.scamp_solve(imageFile,catFile,
		                      file_map.getScampRefCat(fieldName),
		                      filt='r',**kwargs)
		if not keepwcscat:
			os.unlink(catFile)

def make_catalogs(file_map,inputType='sky',**kwargs):
	files = file_map.getFiles(imType='object')
	for imFile in files:
		imageFile = file_map(inputType)(imFile)
		psfFile = file_map('psf',output=True)(imFile)
		if not os.path.exists(psfFile):
			catFile = file_map('wcscat',output=True)(imFile)
			bokphot.sextract(imageFile,catFile,full=False,**kwargs)
			bokphot.run_psfex(catFile,psfFile,**kwargs)
		catFile = file_map('cat',output=True)(imFile)
		bokphot.sextract(imageFile,catFile,psfFile,full=True,**kwargs)

def load_darksky_frames(filt):
	darkSkyFrames = np.loadtxt(os.path.join('config',
	                                        'bokrm_darksky_%s.txt'%filt),
	                           dtype=[('utDate','S8'),('fileName','S35'),
	                                  ('skyVal','f4')])
	# a quick pruning of the repeat images
	return darkSkyFrames[::2]

def create_file_map(obsDb,rawDir,procDir,utds,bands,newfiles,
                    darkskyframes=False,tmpdirin=False,tmpdirout=False):
	# set default file paths
	if rawDir is None:
		rawDir = os.environ['BOK90PRIMERAWDIR']
	if procDir is None:
		procDir = os.path.join(os.environ['BOK90PRIMEOUTDIR'],pipeVersion)
	# create the file manager object
	if newfiles:
		fileMap = ProcessToNewFiles(obsDb,rawDir,procDir) 
	else:
		fileMap = ProcessInPlace(obsDb,rawDir,procDir) 
	if tmpdirin:
		fileMap.setTmpInput()
	if tmpdirout:
		fileMap.setTmpOutput()
	# restrict processing to specified UT dates and filters
	if utds is not None:
		fileMap.setUtDates(utds)
	if bands is not None:
		fileMap.setFilters(bands)
	if darkskyframes:
		# must select a band
		if bands is None or bands not in ['g','i']:
			raise ValueError("Must select a band for dark sky frames (-b)")
		frames = load_darksky_frames(bands)
		fileMap.setFrameList(frames['utDate'],frames['fileName'])
	# create output directories for processed files
	for d in [procDir,fileMap.getCalDir(),fileMap.getDiagDir()]:
		if not os.path.exists(d):
			os.mkdir(d)
	for utd in fileMap.getUtDates():
		utdir = os.path.join(fileMap.procDir,'ut'+utd)
		if not os.path.exists(utdir): os.mkdir(utdir)
		if tmpdirout:
			utdir = os.path.join(fileMap._tmpDir,'ut'+utd)
			if not os.path.exists(utdir): os.mkdir(utdir)
	fileMap.setScampRefCatDir(os.path.join(os.environ['BOK90PRIMEOUTDIR'],
	                                       'scamp_refs'))
	return fileMap

def rmpipe(fileMap,**kwargs):
	redo = kwargs.get('redo',False)
	steps = kwargs.get('steps')
	verbose = kwargs.get('verbose',0)
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
		if kwargs.get('nobiascorr',False):
			biasMap = None
		elif biasMap is None:
			biasMap = get_bias_map(fileMap)
		if kwargs.get('noflatcorr',False):
			flatMap = None
		else:
			flatMap = get_flat_map(fileMap)
		process_all(fileMap,biasMap,flatMap,
		            fixpix=fixpix,
		            norampcorr=kwargs.get('norampcorr'),
		            nocombine=kwargs.get('nocombine'),
		            gain_multiply=not kwargs.get('nogainmul',False),
		            nosavegain=kwargs.get('nosavegain'),
		            prockey=kwargs.get('prockey','CCDPROC'),
		            **pipekwargs)
		timerLog('ccdproc')
	if 'skyflat' in steps:
		make_supersky_flats(fileMap,**pipekwargs)
		timerLog('supersky flats')
	if 'proc2' in steps:
		skyArgs = { k.lstrip('sky'):kwargs[k] 
		                 for k in ['skymethod','skyorder']}
		process_all2(fileMap,skyArgs,
		             noillumcorr=kwargs.get('noillumcorr'),
		             nodarkskycorr=kwargs.get('nodarkskycorr'),
		             prockey=kwargs.get('prockey','CCDPRO2'),
		             save_sky=kwargs.get('savesky'),
		             **pipekwargs)
		timerLog('process2')
	if 'wcs' in steps:
		set_wcs(fileMap,**pipekwargs)
		timerLog('wcs')
	if 'cat' in steps:
		make_catalogs(fileMap,**pipekwargs)
		timerLog('catalog')
	timerLog.dump()

def rmpipe_poormp(fileMap,**kwargs):
	nProc = kwargs.get('processes',1)
	def chunks(l, n):
		nstep = int(round(len(l)/float(n)))
		for i in xrange(0, len(l), nstep):
			yield l[i:i+nstep]
	utdSets = chunks(fileMap.getUtDates(),nProc)
	jobs = []
	for i,utds in enumerate(utdSets):
		fmap = copy(fileMap)
		fmap.setUtDates(utds)
		p = multiprocessing.Process(target=rmpipe,
		                            args=(fmap,),kwargs=kwargs)
		jobs.append(p)
		p.start()

def make_images(file_map,imtype='comb',msktype=None):
	import matplotlib.pyplot as plt
	files = file_map.getFiles(imType='object')
	_fmap = file_map(imtype,False)
	if msktype=='badpix':
		msktype = 'MasterBadPixMask4'
	if msktype==None:
		maskmap = lambda f: None
	else:
		maskmap = file_map(msktype,False)
	imdir = os.path.join(file_map.getProcDir(),'images')
	if not os.path.exists(imdir):
		os.mkdir(imdir)
	plt.ioff()
	for ff in files:
		f = _fmap(ff)
		print ff
		print f
		if not os.path.exists(f):
			continue
		imgfile = os.path.basename(f).replace('.fits','.png')
		imgfile = os.path.join(imdir,imgfile)
		if os.path.exists(imgfile):
			continue
		bokmkimage.make_fov_image_fromfile(f,imgfile,mask=maskmap(ff))
	plt.ion()

if __name__=='__main__':
	import sys
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--band',type=str,default=None,
	                help='band to process (g or i) [default=both]')
	parser.add_argument('-c','--caldir',type=str,default=None,
	                help='set calibration directory')
	parser.add_argument('-f','--file',type=str,default=None,
	                help='file to process [default=all]')
	parser.add_argument('--frames',type=str,default=None,
	                help='frames to process (i1,i2) [default=all]')
	parser.add_argument('-n','--newfiles',action='store_true',
	                help='process to new files (not in-place)')
	parser.add_argument('--obsdb',type=str,default=None,
	                help='location of observations db')
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
	parser.add_argument('-t','--imtype',type=str,default=None,
	                help='specify image type to process')
	parser.add_argument('-u','--utdate',type=str,default=None,
	                help='UT date(s) to process [default=all]')
	parser.add_argument('-v','--verbose',action='count',
	                help='increase output verbosity')
	parser.add_argument('--calccdims',action='store_true',
	                help='generate CCD-combined images for calibration data')
	parser.add_argument('--nobiascorr',action='store_true',
	                help='do not apply bias correction')
	parser.add_argument('--noflatcorr',action='store_true',
	                help='do not apply flat correction')
	parser.add_argument('--norampcorr',action='store_true',
	                help='do not correct bias ramp')
	parser.add_argument('--noillumcorr',action='store_true',
	                help='do not apply illumination correction')
	parser.add_argument('--nodarkskycorr',action='store_true',
	                help='do not apply dark sky flat correction')
	parser.add_argument('--nogainmul',action='store_true',
	                help='do not multiply gain when processing')
	parser.add_argument('--nocombine',action='store_true',
	                help='do not combine into CCD images')
	parser.add_argument('--nosavegain',action='store_true',
	                help='do not save per-image gain balance factors')
	parser.add_argument('--nousepixflat',action='store_true',
	                help='do not use normalized pixel flat')
	parser.add_argument('--noskysub',action='store_true',
	                help='do not perform sky subtraction')
	parser.add_argument('--skymethod',type=str,default='polynomial',
	                help='sky subtraction method ([polynomial]|spline)')
	parser.add_argument('--skyorder',type=int,default=1,
	                help='sky subtraction order [default: 1 (linear)]')
	parser.add_argument('--savesky',action='store_true',
	                help='save sky background fit')
	parser.add_argument('--darkskyframes',action='store_true',
	                help='load only the dark sky frames')
	parser.add_argument('--prockey',type=str,default=None,
	                help='set new header key for ccdproc')
	parser.add_argument('--tmpdirin',action='store_true',
	                help='read files from temporary directory')
	parser.add_argument('--tmpdirout',action='store_true',
	                help='write files to temporary directory')
	parser.add_argument('--images',type=str,default=None,
	                help='make png images (imtype,[msktype]) '
	                     'instead of processing')
	parser.add_argument('--makerampcorr',action='store_true',
	                help='make ramp correction image instead of processing')
	parser.add_argument('--makeillumcorr',action='store_true',
	                help='make illumination correction image '
	                     'instead of processing')
	parser.add_argument('--wcscheck',action='store_true',
	                help='make astrometry diagnostic files')
	args = parser.parse_args()
	if args.obsdb is None:
		obsDb = Table.read(os.path.join('config','sdssrm-bok2014.fits'))
	else:
		obsDb = Table.read(args.obsdb)
	if args.utdate is None:
		utds = None
	else:
		utds = args.utdate.split(',')
	if args.steps is None:
		if args.stepto is None:
			steps = all_process_steps
		else:
			steps = all_process_steps[:all_process_steps.index(args.stepto)+1]
	elif args.steps == 'proc':
		steps = ['oscan','proc1','proc2']
	else:
		steps = args.steps.split(',')
	verbose = 0 if args.verbose is None else args.verbose
	# set up the data map
	fileMap = create_file_map(obsDb,args.rawdir,args.output,
	                          utds,args.band,args.newfiles,args.darkskyframes,
	                          args.tmpdirin,args.tmpdirout)
	if args.frames is not None:
		fileMap.setFrames(tuple([int(_f) for _f in args.frames.split(',')]))
	elif args.file is not None:
		fileMap.setFile(args.file)
	if args.imtype is not None:
		fileMap.setImageType(args.imtype)
	if args.caldir is not None:
		fileMap.setCalDir(os.path.join(args.caldir,'cals'))
		fileMap.setDiagDir(os.path.join(args.caldir,'diagnostics'))
	# convert command-line arguments into dictionary
	opts = vars(args)
	kwargs = { k : opts[k] for k in opts if opts[k] != None }
	kwargs['steps'] = steps
	# run pipeline processes
	if args.images is not None:
		make_images(fileMap,*args.images.split(','))
	elif args.makerampcorr:
		bokrampcorr.make_rampcorr_image(fileMap)#,**kwargs)
	elif args.makeillumcorr:
		bokillumcorr.make_illumcorr_image(fileMap)#,**kwargs)
	elif args.wcscheck:
		files = [fileMap('sky')(f) for f in fileMap.getFiles('object')]
		bokgnostic.run_scamp_diag(files)
	elif args.processes > 1:
		rmpipe_poormp(fileMap,**kwargs)
	else:
		rmpipe(fileMap,**kwargs)

