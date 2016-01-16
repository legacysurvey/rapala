#!/usr/bin/env python

import os
import re
import glob
import shutil
from copy import copy
import multiprocessing
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio
from astropy.table import Table

from .bokoscan import BokOverscanSubtract
from .badpixels import build_mask_from_flat
from . import bokio
from . import bokutil
from . import bokproc
from . import bokrampcorr
from . import bokillumcorr
from . import bokastrom
from . import bokphot
from . import bokgnostic

all_process_steps = ['oscan','bias2d','flat2d','bpmask',
                     'proc1','comb','skyflat','proc2','sky','wcs','cat']

default_filenames = {
  'oscan':'','bias':'_b','proc1':'_p','comb':'_c',
  'pass1cat':'.cat1','skymask':'.skymsk','skyfit':'.sky',
  'sky':'_s','proc2':'_q','wcscat':'.wcscat',
  'cat':'.cat','psf':'.psf'
}

class SimpleFileNameMap(bokio.FileNameMap):
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
				if not os.path.exists(fpath):
					# or fpacked
					fpath = os.path.join(self.rawDir,fn+'.fz')
			return fpath
		else:
			return os.path.join(self.procDir,fn)

class BokDataManager(object):
	def __init__(self,obsDb,rawDir,procDir):
		self.rawDir = rawDir
		self.procDir = procDir
		self.calDir = os.path.join(self.procDir,'cals')
		self.diagDir = os.path.join(self.procDir,'diagnostics')
		self.master = {}
		self.obsDb = obsDb
		self.allUtDates = np.unique(self.obsDb['utDate'])
		self.utDates = self.allUtDates
		self.allFilt = np.unique(obsDb['filter'][obsDb['imType']=='object'])
		self.filt = self.allFilt
		self.frames = None
		self.frameList = None
		self.imType = None
		self._curUtDate = None
		self._curFilt = None
		self._tmpInput = False
		self._tmpOutput = False
		self._tmpDir = os.path.join(self.procDir,'tmp')
		self.refCatDir = None
		self.fileSuffixes = default_filenames
		self.setInPlace(True)
		self.firstStep = None
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
		for utdir in self.getUtDirs():
			utdir = os.path.join(self._tmpDir,utdir)
			if not os.path.exists(utdir): os.mkdir(utdir)
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
	def getUtDirs(self):
		return sorted(np.unique(self.obsDb['utDir']))
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
	def setFileList(self,utDates,fileNames):
		frames = [ np.where( (self.obsDb['utDate']==utd) &
		                     (self.obsDb['fileName']==f) )[0][0]
		             for utd,f in zip(utDates,fileNames) ]
		self.setFrameList(frames)
	def setFrameList(self,frames):
		self.frameList = np.array(frames)
	def setImageType(self,imType):
		self.imType = imType
	def setScampRefCatDir(self,refCatDir):
		self.refCatDir = refCatDir
		if not os.path.exists(self.refCatDir):
			os.mkdir(self.refCatDir)
	def getScampRefCat(self,fieldName):
		if self.refCatDir is None:
			return None
		refCatFn = 'scampref_%s.cat' % fieldName
		return os.path.join(self.refCatDir,refCatFn)
	def setMaster(self,calType,fits=None,byFilter=False):
		if byFilter:
			self.master[calType] = {}
			for filt in self.filt:
				self.master[calType][filt] = {'fits':None,
				                     'fileName':'%s_%s.fits'%(calType,filt)}
		else:
			self.master[calType] = {'fits':fits,'fileName':calType+'.fits'}
	def getMaster(self,calType,filt=None,name=False,load=True):
		master = self.master[calType]
		if filt is not None:
			master = master[filt]
		if name:
			return master['fileName']
		else:
			if load and master['fits'] is None:
				master['fits'] = fitsio.FITS(os.path.join(self.calDir,
				                                      master['fileName']))
			return master['fits']
	def setInPlace(self,inPlace):
		self.inPlace = inPlace
		if self.inPlace:
			self.filesToMap = ['oscan','pass1cat','skymask','skyfit',
			                   'wcscat','cat','psf']
		else:
			self.filesToMap = self.fileSuffixes.keys()
	def __call__(self,t):
		procDir = self.procDir
		if t in all_process_steps:
			if self.firstStep is None:
				self.firstStep = t
				if self._tmpInput:
					procDir = self._tmpDir
			elif t != self.firstStep:
				if self._tmpOutput:
					procDir = self._tmpDir
			else:
				if self._tmpInput:
					procDir = self._tmpDir
		if t.startswith('Master'):
			_t = t.lstrip('Master')
			if self.master[_t]['fits'] is None:
				fn = os.path.join(self.calDir,self.master[_t]['fileName'])
				self.master[_t]['fits'] = fitsio.FITS(fn)
			return self.master[_t]['fits']
		elif t == 'raw':
			# special case to handle moving files from raw directory
			return SimpleFileNameMap(self.rawDir,procDir,fromRaw=True)
		elif t in self.filesToMap:
			return SimpleFileNameMap(self.rawDir,procDir,
			                         self.fileSuffixes[t])
		else:
			return SimpleFileNameMap(self.rawDir,procDir)
	def getFiles(self,imType=None,utd=None,filt=None,
	             im_range=None,exclude_objs=None,
	             with_objnames=False,with_frames=False,as_sequences=False):
		try:
			# if the observations database has a good flag use it
			good_ims = self.obsDb['good']
		except:
			good_ims = np.ones(len(self.obsDb),dtype=bool)
		file_sel = good_ims.copy()
		# select on UT date(s)
		if utd is not None:
			utds = [utd] if type(utd) is str else utd
		elif self._curUtDate is not None:
			utds = [self._curUtDate]
		else:
			utds = self.getUtDates() 
		if not np.array_equal(utds,self.allUtDates):
			isUtd = np.zeros_like(file_sel)
			for utd in utds:
				isUtd |= ( self.obsDb['utDate'] == utd )
			file_sel &= isUtd
		# restrict on image type
		if imType is not None:
			file_sel &= ( self.obsDb['imType'] == imType )
		elif self.imType is not None:
			file_sel &= ( self.obsDb['imType'] == self.imType )
		# restrict on filter
		if filt is not None:
			f = filt
		elif self._curFilt is not None:
			f = self._curFilt
		else:
			f = self.filt
		if not np.all(np.in1d(self.allFilt,f)):
			isFilt = np.in1d(self.obsDb['filter'],f)
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
			if as_sequences:
				splits = np.where(np.diff(ii)>1)[0]
				if len(splits)==0:
					return [ files ]
				else:
					# look for instances where the sequence is broken simply
					# because of a few bad images
					_splits = []
					for s in splits:
						if not np.all(~good_ims[ii[s]+1:ii[s+1]]):
							_splits.append(s)
					if len(_splits)==0:
						return [ files ]
					else:
						return np.split(files,np.array(_splits)+1)
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

def get_bias_map(dataMap):
	biasMap = {}
	biasPattern = os.path.join(dataMap.getCalDir(),'bias_*.fits')
	biasFiles = sorted(glob.glob(biasPattern))
	bias2utd = ".*bias_(\d+).*"
	biasUtds = np.array([int(re.match(bias2utd,fn).groups()[0])
	                       for fn in biasFiles])
	for utd in dataMap.iterUtDates():
		j = np.argmin(np.abs(int(utd)-biasUtds))
		biasFile = biasFiles[j]
		files = dataMap.getFiles()
		if files is None:
			continue
		for f in files:
			biasMap[f] = biasFile
	return biasMap

def get_flat_map(dataMap):
	flatMap = {}
	flatPattern = os.path.join(dataMap.getCalDir(),'flat_????????_*.fits')
	flatFiles = sorted(glob.glob(flatPattern))
	flat2utdfilt = ".*flat_(\d+)_(\w+)_.*"
	utdfilt = [ re.match(flat2utdfilt,fn).groups() for fn in flatFiles]
	flatUtds = np.array([int(utd) for utd,filt in utdfilt])
	flatFilt = np.array([filt for utd,filt in utdfilt])
	for filt in dataMap.iterFilters():
		for utd in dataMap.iterUtDates():
			files = dataMap.getFiles()
			if files is None:
				continue
			jj = np.where(flatFilt==filt)[0]
			try:
				j = np.argmin(np.abs(int(utd)-flatUtds[jj]))
			except ValueError:
				raise ValueError('no match for %s:%s in %s:%s' % 
				                 (utd,filt,flatUtds,flatFilt))
			flatFile = flatFiles[jj[j]]
			for f in files:
				flatMap[f] = flatFile
	return flatMap

def makeccd4image(dataMap,inputFile,outputFile=None,**kwargs):
	if outputFile is None:
		ccd4map = bokio.FileNameMap(dataMap.getCalDir(),'_4ccd')
	else:
		ccd4map = lambda f: outputFile
	bokproc.combine_ccds([inputFile,],output_map=ccd4map,**kwargs)

def overscan_subtract(dataMap,**kwargs):
	oscanSubtract = BokOverscanSubtract(input_map=dataMap('raw'),
                                        output_map=dataMap('oscan'),
                                        **kwargs)
	oscanSubtract.process_files(dataMap.getFiles())

def make_2d_biases(dataMap,nSkip=2,reject='sigma_clip',
                   writeccdim=False,**kwargs):
	biasStack = bokproc.BokBiasStack(input_map=dataMap('oscan'),
	                                 reject=reject,
                                     **kwargs)
	for utd in dataMap.iterUtDates():
		bias_seqs = dataMap.getFiles(imType='zero',as_sequences=True)
		if bias_seqs is None:
			continue
		for biasNum,biasFiles in enumerate(bias_seqs,start=1):
			if len(biasFiles) < 5: # XXX hardcoded
				continue
			biasFile = os.path.join(dataMap.getCalDir(),
			                        'bias_%s_%d.fits' % (utd,biasNum))
			biasStack.stack(biasFiles[nSkip:],biasFile)
			if kwargs.get('verbose',0) >= 1:
				print '2DBIAS: ',biasFile
				print '\n'.join(biasFiles[nSkip:])
			if writeccdim:
				makeccd4image(dataMap,biasFile,**kwargs)

def make_dome_flats(dataMap,bias_map,
                    nSkip=1,reject='sigma_clip',writeccdim=False,
	                usepixflat=True,debug=False,**kwargs):
	bias2Dsub = bokproc.BokCCDProcess(input_map=dataMap('oscan'),
	                                  output_map=dataMap('bias'),
	                                  bias=bias_map,
	                                  header_key='BIAS2D',
	                                  **kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=dataMap('bias'),
	                                     **kwargs)
	if usepixflat:
		if debug:
			ffmap = SimpleFileNameMap(None,dataMap.procDir,'_fit')
			bfmap = SimpleFileNameMap(None,dataMap.procDir,'_binned')
		else:
			ffmap,bfmap = None,None
		normFlat = bokproc.NormalizeFlat(_normed_flat_fit_map=ffmap,
		                                 _binned_flat_map=bfmap,**kwargs)
	for utd in dataMap.iterUtDates():
		for filt in dataMap.iterFilters():
			flat_seqs = dataMap.getFiles(imType='flat',as_sequences=True)
			if flat_seqs is None:
				continue
			for flatNum,flatFiles in enumerate(flat_seqs,start=1):
				if len(flatFiles) < 5: # XXX hardcoded
					continue
				flatFile = os.path.join(dataMap.getCalDir(),
			                        'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
				bias2Dsub.process_files(flatFiles)
				flatStack.stack(flatFiles[nSkip:],flatFile)
				if usepixflat:
					if debug:
						shutil.copy(flatFile,
						            flatFile.replace('.fits','_raw.fits'))
					normFlat.process_files([flatFile])
				if kwargs.get('verbose',0) >= 1:
					print 'DOMEFLAT: ',flatFile
					print '\n'.join(flatFiles[nSkip:])
				if writeccdim:
					makeccd4image(dataMap,flatFile,**kwargs)

def make_bad_pixel_masks(dataMap,**kwargs):
	for utd in dataMap.iterUtDates():
		for filt in dataMap.iterFilters():
			flatNum = 1
			flatFile = os.path.join(dataMap.getCalDir(),
			                      'flat_%s_%s_%d.fits' % (utd,filt,flatNum))
			if not os.path.exists(flatFile):
				continue
			hdr = fitsio.read_header(flatFile,0)
			if 'NORMFLT' not in hdr:
				# need to normalize the flat before making mask
				_map = SimpleFileNameMap(None,dataMap.procDir,'_normed')
				ffmap = SimpleFileNameMap(None,dataMap.procDir,'_fit')
				bfmap = SimpleFileNameMap(None,dataMap.procDir,'_binned')
				normFlat = bokproc.NormalizeFlat(output_map=_map,
				                            _normed_flat_fit_map=ffmap,
				                            _binned_flat_map=bfmap,**kwargs)
				normFlat.process_files([flatFile])
				flatFile = _map(flatFile)
			bpMaskFile = os.path.join(dataMap.getCalDir(),
			                       dataMap.getMaster('BadPixMask',name=True))
			bpMask4File = os.path.join(dataMap.getCalDir(),
			                       dataMap.getMaster('BadPixMask4',name=True))
			if kwargs.get('verbose',0) >= 1:
				print 'generating bad pixel mask ',bpMaskFile,
				print ' from ',flatFile
			build_mask_from_flat(flatFile,bpMaskFile)#,**kwargs)
			makeccd4image(dataMap,bpMaskFile,outputFile=bpMask4File,**kwargs)
			break

def balance_gains(dataMap,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=dataMap('proc1'),
	                                     mask_map=dataMap('MasterBadPixMask'),
	                                                **kwargs)
	gainMap = {'corrections':{},'skyvals':{}}
	for utd in dataMap.iterUtDates():
		for filt in dataMap.iterFilters():
			files = dataMap.getFiles(imType='object')
			if files is None:
				continue
			diagfile = os.path.join(dataMap.getDiagDir(),
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

def process_all(dataMap,bias_map,flat_map,
                fixpix=False,norampcorr=False,
                nocombine=False,prockey='CCDPROC',**kwargs):
	# 1. basic processing (bias and flat-field correction, fixpix, 
	#    nominal gain correction
	ramp = None if norampcorr else dataMap('MasterBiasRamp')
	proc = bokproc.BokCCDProcess(input_map=dataMap('oscan'),
	                             output_map=dataMap('proc1'),
	                             mask_map=dataMap('MasterBadPixMask'),
	                             header_key=prockey,
	                             bias=bias_map,flat=flat_map,
	                             ramp=ramp,illum=None,darksky=None,
	                             fixpix=fixpix,**kwargs)
	files = dataMap.getFiles(imType='object')
	if files is None:
		return
	proc.process_files(files)
	if nocombine:
		return
	# 2. balance gains using background counts
	gainMap = balance_gains(dataMap,**kwargs)
	# 3. combine per-amp images (16) into CCD images (4)
	bokproc.combine_ccds(files,
	                     input_map=dataMap('proc1'), 
	                     output_map=dataMap('comb'),
	                     gain_map=gainMap,
	                     **kwargs)

def make_supersky_flats(dataMap,**kwargs):
	caldir = dataMap.getCalDir()
	stackin = dataMap('sky') # XXX
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=dataMap('skymask'),
	                    exposure_time_map=bokio.FileNameMap(caldir,'.exp'),
	                       raw_stack_file=bokio.FileNameMap(caldir,'_raw'),
	                                        header_bad_key='BADSKY')
	skyFlatStack.set_badpixelmask(dataMap('MasterBadPixMask4'))
	for filt in dataMap.iterFilters():
		if True:
			files = dataMap.getFiles(imType='object')
			if files is not None:
				fn = dataMap.getMaster('SkyFlat',filt=filt,name=True)
				outfn = os.path.join(dataMap.getCalDir(),fn)
				skyFlatStack.stack(files,outfn)
		if False:
			# this does a sky flat each night
			for utd in dataMap.iterUtDates():
				files = dataMap.getFiles(imType='object')
				if files is None:
					continue
				outfn = os.path.join(dataMap.getCalDir(),
				                     'skyflat_%s_%s.fits' % (utd,filt))
				skyFlatStack.stack(files,outfn)

def process_all2(dataMap,skyArgs,noillumcorr=False,nodarkskycorr=False,
                 noskysub=False,prockey='CCDPRO2',save_sky=False,
                 **kwargs):
	#
	# Second round flat-field corrections
	#
	files,ii = dataMap.getFiles(imType='object',with_frames=True)
	if files is None:
		return
	if noillumcorr:
		illum_map = None
	else:
		illum_map = {}
		for i,f in zip(ii,files):
			illum_map[f] = dataMap.getMaster('Illumination',
			                                 dataMap.obsDb['filter'][i])
	if nodarkskycorr:
		darksky_map = None
	else:
		darksky_map[f] = {}
		for i,f in zip(ii,files):
			darksky_map[f] = dataMap.getMaster('SkyFlat',
			                                   dataMap.obsDb['filter'][i])
	#
	proc = bokproc.BokCCDProcess(input_map=dataMap('comb'),
	                             output_map=dataMap('proc2'),
	                             mask_map=dataMap('MasterBadPixMask4'),
	                             header_key=prockey,
	                             gain_multiply=False,bias=None,flat=None,
	                             ramp=None,fixpix=False,
	                             illum=illum_map,darksky=darksky_map,
	                             **kwargs)
	proc.process_files(files)
	if noskysub:
		return
	#
	# Sky subtraction
	#
	# Generate sky masks by agressively masking objects
	skyFlatMask = bokproc.BokGenerateSkyFlatMasks(
	                                    input_map=dataMap('proc2'),
	                                    output_map=dataMap('skymask'),
	                                    mask_map=dataMap('MasterBadPixMask4'))
	files = dataMap.getFiles(imType='object')
	skyFlatMask.process_files(files)
	skyfitmap = dataMap('skyfit') if save_sky else None
	skySub = bokproc.BokSkySubtract(input_map=dataMap('proc2'),
	                                output_map=dataMap('sky'),
	                                mask_map=dataMap('skymask'),
	                                skyfit_map=skyfitmap,**skyArgs)
	skySub.add_mask(dataMap('MasterBadPixMask4'))
	skySub.process_files(files)

def set_wcs(dataMap,inputType='sky',keepwcscat=True,**kwargs):
	filesAndFields = dataMap.getFiles(imType='object',with_objnames=True)
	for imFile,fieldName in zip(*filesAndFields):
		imageFile = dataMap(inputType)(imFile)
		catFile = dataMap('wcscat')(imFile)
		bokphot.sextract(imageFile,catFile,**kwargs)
		bokastrom.scamp_solve(imageFile,catFile,
		                      dataMap.getScampRefCat(fieldName),
		                      filt='r',**kwargs)
		if not keepwcscat:
			os.unlink(catFile)

def make_catalogs(dataMap,inputType='sky',**kwargs):
	files = dataMap.getFiles(imType='object')
	for imFile in files:
		imageFile = dataMap(inputType)(imFile)
		psfFile = dataMap('psf')(imFile)
		if not os.path.exists(psfFile):
			catFile = dataMap('wcscat')(imFile)
			bokphot.sextract(imageFile,catFile,full=False,**kwargs)
			bokphot.run_psfex(catFile,psfFile,**kwargs)
		catFile = dataMap('cat')(imFile)
		bokphot.sextract(imageFile,catFile,psfFile,full=True,**kwargs)

def rmpipe(dataMap,**kwargs):
	redo = kwargs.get('redo',False)
	steps = kwargs.get('steps')
	debug = kwargs.get('debug',False)
	verbose = kwargs.get('verbose',0)
	pipekwargs = {'clobber':redo,'verbose':verbose}
	# fixpix is sticking nan's into the images in unmasked pixels (???)
	fixpix = False #True
	writeccdims = kwargs.get('calccdims',False)
	timerLog = bokutil.TimerLog()
	biasMap = None
	if 'oscan' in steps:
		overscan_subtract(dataMap,**pipekwargs)
		timerLog('overscans')
	if 'bias2d' in steps:
		make_2d_biases(dataMap,writeccdim=writeccdims,**pipekwargs)
		timerLog('2d biases')
	if 'flat2d' in steps:
		biasMap = get_bias_map(dataMap)
		if debug:
			pass
		make_dome_flats(dataMap,biasMap,writeccdim=writeccdims,
		                usepixflat=not kwargs.get('nousepixflat',False),
		                debug=debug,**pipekwargs)
		timerLog('dome flats')
	if 'bpmask' in steps:
		make_bad_pixel_masks(dataMap)
		timerLog('bad pixel masks')
	if 'proc1' in steps or 'comb' in steps:
		if kwargs.get('nobiascorr',False):
			biasMap = None
		elif biasMap is None:
			biasMap = get_bias_map(dataMap)
		if kwargs.get('noflatcorr',False):
			flatMap = None
		else:
			flatMap = get_flat_map(dataMap)
		process_all(dataMap,biasMap,flatMap,
		            fixpix=fixpix,
		            norampcorr=kwargs.get('norampcorr'),
		            nocombine=kwargs.get('nocombine'),
		            gain_multiply=not kwargs.get('nogainmul',False),
		            nosavegain=kwargs.get('nosavegain'),
		            prockey=kwargs.get('prockey','CCDPROC'),
		            **pipekwargs)
		timerLog('ccdproc')
	if 'skyflat' in steps:
		make_supersky_flats(dataMap,**pipekwargs)
		timerLog('supersky flats')
	if 'proc2' in steps:
		skyArgs = { k.lstrip('sky'):kwargs[k] 
		                 for k in ['skymethod','skyorder']}
		process_all2(dataMap,skyArgs,
		             noillumcorr=kwargs.get('noillumcorr'),
		             nodarkskycorr=kwargs.get('nodarkskycorr'),
		             prockey=kwargs.get('prockey','CCDPRO2'),
		             save_sky=kwargs.get('savesky'),
		             **pipekwargs)
		timerLog('process2')
	if 'wcs' in steps:
		set_wcs(dataMap,**pipekwargs)
		timerLog('wcs')
	if 'cat' in steps:
		make_catalogs(dataMap,**pipekwargs)
		timerLog('catalog')
	timerLog.dump()

def rmpipe_poormp(dataMap,**kwargs):
	nProc = kwargs.get('processes',1)
	def chunks(l, n):
		nstep = int(round(len(l)/float(n)))
		for i in xrange(0, len(l), nstep):
			yield l[i:i+nstep]
	utdSets = chunks(dataMap.getUtDates(),nProc)
	jobs = []
	for i,utds in enumerate(utdSets):
		fmap = copy(dataMap)
		fmap.setUtDates(utds)
		p = multiprocessing.Process(target=rmpipe,
		                            args=(fmap,),kwargs=kwargs)
		jobs.append(p)
		p.start()

def make_images(dataMap,imtype='comb',msktype=None):
	import matplotlib.pyplot as plt
	files = dataMap.getFiles(imType='object')
	_fmap = dataMap(imtype,False)
	if msktype=='badpix':
		msktype = 'MasterBadPixMask4'
	if msktype==None:
		maskmap = lambda f: None
	else:
		maskmap = dataMap(msktype,False)
	imdir = os.path.join(dataMap.getProcDir(),'images')
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

def init_file_args(parser):
	parser.add_argument('-b','--band',type=str,default=None,
	                help='band to process (g or i) [default=both]')
	parser.add_argument('-c','--caldir',type=str,default=None,
	                help='set calibration directory')
	parser.add_argument('-f','--file',type=str,default=None,
	                help='file to process [default=all]')
	parser.add_argument('--frames',type=str,default=None,
	                help='frames to process (i1,i2) [default=all]')
	parser.add_argument('--obsdb',type=str,default=None,
	                help='location of observations db')
	parser.add_argument('-o','--output',type=str,default=None,
	                help='output directory [default=$BOK90PRIMEOUTDIR]')
	parser.add_argument('-r','--rawdir',type=str,default=None,
	                help='raw data directory [default=$BOK90PRIMERAWDIR]')
	parser.add_argument('-R','--redo',action='store_true',
	                help='redo (overwrite existing files)')
	parser.add_argument('-t','--imtype',type=str,default=None,
	                help='specify image type to process')
	parser.add_argument('-u','--utdate',type=str,default=None,
	                help='UT date(s) to process [default=all]')
	return parser

def _load_obsdb(obsdb):
	obsDb = Table.read(obsdb)
	# found that when topcat writes FITS tables it adds whitespace to str
	# columns, strip them here (makes a copy but oh well)
	for k in ['utDate','utDir','fileName','imType','filter','objName']:
		obsDb[k] = np.char.rstrip(obsDb[k])
	return obsDb

def init_data_map(args,create_dirs=True):
	#
	obsDb = _load_obsdb(args.obsdb)
	# set up the data map
	dataMap = BokDataManager(obsDb,args.rawdir,args.output)
	#
	if args.utdate is not None:
		dataMap.setUtDates(args.utdate.split(','))
	#
	if args.band is not None:
		dataMap.setFilters(args.band.split(','))
	#
	if args.caldir is not None:
		dataMap.setCalDir(os.path.join(args.caldir,'cals'))
		dataMap.setDiagDir(os.path.join(args.caldir,'diagnostics'))
	#
	if args.frames is not None:
		_frames = []
		for _f in args.frames.split(','):
			if '-' in _f:
				i1,i2 = _f.split('-')
				_frames.extend(range(int(i1),int(i2)+1))
			else:
				_frames.append(int(_f))
		dataMap.setFrameList(_frames)
	elif args.file is not None:
		dataMap.setFile(args.file)
	if args.imtype is not None:
		dataMap.setImageType(args.imtype)
	#
	if create_dirs:
		for d in [dataMap.procDir,dataMap.getCalDir(),dataMap.getDiagDir()]:
			if not os.path.exists(d):
				os.mkdir(d)
		for _utdir in dataMap.getUtDirs():
			utdir = os.path.join(dataMap.procDir,_utdir)
			if not os.path.exists(utdir): os.mkdir(utdir)
	return dataMap

def set_master_cals(dataMap):
	dataMap.setMaster('BadPixMask')
	dataMap.setMaster('BadPixMask4')
	dataMap.setMaster('BiasRamp')
	dataMap.setMaster('Illumination',byFilter=True)
	dataMap.setMaster('SkyFlat',byFilter=True)
	return dataMap

def init_pipeline_args(parser):
	parser.add_argument('--debug',action='store_true',
	                help='save additional debugging files')
	parser.add_argument('-n','--newfiles',action='store_true',
	                help='process to new files (not in-place)')
	parser.add_argument('-p','--processes',type=int,default=1,
	                help='number of processes to use [default=single]')
	parser.add_argument('-s','--steps',type=str,default=None,
	                help='processing steps to execute [default=all]')
	parser.add_argument('-S','--stepto',type=str,default=None,
	                help='process until this step [default=last]')
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
	return parser

def run_pipe(dataMap,args):
	if args.newfiles:
		dataMap.setInPlace(False)
	if args.tmpdirin:
		dataMap.setTmpInput()
	if args.tmpdirout:
		dataMap.setTmpOutput()
	if args.steps is None:
		if args.stepto is None:
			steps = all_process_steps
		else:
			steps = all_process_steps[:all_process_steps.index(args.stepto)+1]
	elif args.steps == 'proc':
		steps = ['oscan','proc1','proc2']
	elif args.steps == 'all':
		steps = ['oscan','proc1','proc2','wcs','cat']
	else:
		steps = args.steps.split(',')
	verbose = 0 if args.verbose is None else args.verbose
	# convert command-line arguments into dictionary
	opts = vars(args)
	kwargs = { k : opts[k] for k in opts if opts[k] != None }
	kwargs['steps'] = steps
	# run pipeline processes
	if args.images is not None:
		make_images(dataMap,*args.images.split(','))
	elif args.makerampcorr:
		bokrampcorr.make_rampcorr_image(dataMap)#,**kwargs)
	elif args.makeillumcorr:
		bokillumcorr.make_illumcorr_image(dataMap)#,**kwargs)
	elif args.wcscheck:
		files = [dataMap('sky')(f) for f in dataMap.getFiles('object')]
		bokgnostic.run_scamp_diag(files)
	elif args.processes > 1:
		rmpipe_poormp(dataMap,**kwargs)
	else:
		rmpipe(dataMap,**kwargs)

