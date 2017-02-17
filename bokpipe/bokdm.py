#!/usr/bin/env python

import os
import pickle
import numpy as np
from numpy.core.defchararray import add as char_add
from astropy.table import Table

from .bokio import FileNameMap,IdentityNameMap
from .bokutil import FakeFITS,array_stats,stats_region,load_mask

##############################################################################
#                                                                            #
# Calibration Database                                                       #
#   identify sequences of calibration frames and store in a database         #
#                                                                            #
##############################################################################

def find_cal_sequences(obsDb,imType,byFilt=True,
                       minLen=5,maxDt=None,filts=None):
	# first group all the data by night
	t = obsDb.group_by('utDir')
	# frameIndex isn't actually a straight table index, but rather a unique
	# id. adding a running index makes the group sorting much clearer.
	t['ii'] = np.arange(len(t))
	calseqs = []
	if byFilt:
		tgroups = ['imType','filter','expTime']
	else:
		tgroups = ['imType','expTime']
	for ut in t.groups:
		if filts is None or not byFilt:
			isfilt = True
		else:
			isfilt = np.in1d(ut['filter'],filts)
		iscal = np.where((ut['imType']==imType) & isfilt)[0]
		if len(iscal)==0:
			continue
		ut_type = ut[iscal].group_by(tgroups)
		for utt in ut_type.groups:
			if len(utt) < minLen:
				continue
			ii = np.arange(len(utt))
			seqs = np.split(ii,1+np.where(np.diff(utt['ii'])>1)[0])
			seqs = [ list(utt['ii'][s]) for s in seqs if len(s) >= minLen ]
			calseqs.extend(seqs)
	return calseqs

def caldb_store(calDbFile,obsDb,imType,seqs,useFilt=True):
	prefixes = {'zero':'Bias','flat':'DomeFlat',
	            'illum':'Illum','fringe':'Fringe','skyflat':'SkyFlat'}
	pfx = prefixes.get(imType,imType)
	try:
		with open(calDbFile,"rb") as caldbf:
			calDb = pickle.load(caldbf)
	except:
		calDb = {}
	rv = []
	calDb.setdefault(imType,[])
	for seq in seqs:
		utd = obsDb['utDate'][seq[0]]
		filt = obsDb['filter'][seq[0]]
		mjd = np.mean(obsDb['mjdStart'][seq])
		utt = obsDb['utObs'][seq[0]][:5].replace(':','')
		fileName = '%s%8s%4s' % (pfx,utd,utt)
		if useFilt:
			fileName += filt
		# filter out any existing entries with same output filename
		# to avoid duplicates
		calDb[imType] = filter(lambda it: it[0] != fileName, calDb[imType])
		calDb[imType].append((fileName,utd,mjd,filt,seq))
		rv.append(fileName)
	with open(calDbFile,"wb") as caldbf:
		pickle.dump(calDb,caldbf)
	return rv

def load_caldb(calDbFile):
	with open(calDbFile,"rb") as caldbf:
		calDb = pickle.load(caldbf)
	return calDb

def write_caldb(calDbFile,calDb):
	with open(calDbFile,"wb") as caldbf:
		pickle.dump(calDb,caldbf)

def init_cal_db(calDbFile,obsDb,filts,overwrite=False):
	if overwrite:
		try:
			os.unlink(calDbFile)
		except:
			pass
	biasSeqs = find_cal_sequences(obsDb,'zero',byFilt=False,maxDt=10)
	flatSeqs = find_cal_sequences(obsDb,'flat',byFilt=True,maxDt=10,
	                              filts=filts)
	caldb_store(calDbFile,obsDb,'zero',biasSeqs,useFilt=False)
	caldb_store(calDbFile,obsDb,'flat',flatSeqs)
	return load_caldb(calDbFile)

##############################################################################
#                                                                            #
# SimpleFileNameMap                                                          #
#   default mapping of input file names to output file names, mapping from   #
#   raw data directory to output directory, and adding suffixes when keeping #
#   images created during intermediate steps.                                #
#                                                                            #
##############################################################################

default_filenames = {
  'oscan':'','bias':'_b','proc1':'_p','comb':'_c','weight':'.wht',
  'imgmask':'.dq','pass1cat':'.cat1','skymask':'.skymsk','skyfit':'.sky',
  'sky':'_s','proc2':'_q','wcscat':'.wcscat',
  'cat':'.cat','psf':'.psf'
}

class SimpleFileNameMap(FileNameMap):
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

#####

class BokCalibrator(object):
	def setTarget(self,f):
		raise NotImplementedError
	def getImage(self,extn):
		raise NotImplementedError
	def getFileName(self):
		raise NotImplementedError
	def __getitem__(self,extn):
		return self.getImage(extn)

class NullCalibrator(BokCalibrator):
	def setTarget(self,f):
		pass
	def getImage(self,extn):
		return None
	def getFileName(self):
		return None

class MasterCalibrator(BokCalibrator):
	def __init__(self,masterFile):
		super(MasterCalibrator,self).__init__()
		self.masterFile = masterFile
		self.masterFits = None
	def _load_fits(self):
		if self.masterFits is None:
			print 'Loading master cal ',os.path.basename(self.masterFile)
			self.masterFits = FakeFITS(self.masterFile)
	def setTarget(self,f):
		pass
	def getImage(self,extn):
		self._load_fits()
		return self.masterFits[extn]
	def getFileName(self):
		return self.masterFile
	def __call__(self,f):
		self._load_fits()
		return self.masterFits

class CalibratorMap(BokCalibrator):
	def __init__(self,obsDb,calTab,nameMap=None,allowMissing=False):
		super(CalibratorMap,self).__init__()
		self.allowMissing = allowMissing
		self.currentFile = None
		self.currentFits = None
		self.calMap = {}
		if nameMap is None:
			nameMap = IdentityNameMap
		self.nameMap = nameMap
		domap = np.in1d(obsDb['imType'],['flat','object'])
		if 'filter' in calTab.colnames:
			for filt in np.unique(calTab['filter']):
				ii = np.where(calTab['filter'] == filt)[0]
				jj = np.where((obsDb['filter'] == filt) & domap)[0]
				for j in jj:
					k = os.path.join(obsDb['utDir'][j],obsDb['fileName'][j])
					dt = np.abs(obsDb['mjdStart'][j]-calTab['mjd'][ii])
					i = np.argmin(dt)
					self.calMap[k] = self.nameMap(calTab['fileName'][ii[i]])
		else:
			jj = np.where(domap)[0]
			for j in jj:
				i = np.argmin(np.abs(obsDb['mjdStart'][j]-calTab['mjd']))
				k = os.path.join(obsDb['utDir'][j],obsDb['fileName'][j])
				self.calMap[k] = self.nameMap(calTab['fileName'][i])
	def setTarget(self,f):
		try:
			cal = self.calMap[f]
		except KeyError:
			if self.allowMissing:
				return
			else:
				raise ValueError("image %s has no calibration" % f)
		if cal != self.currentFile:
			calfn = os.path.basename(cal)
			prevfn = '<None>' if not self.currentFile \
			                      else os.path.basename(self.currentFile) 
			print 'reset cal %s %s %s' % (f,calfn,prevfn)
			if self.currentFits:
				self.currentFits.close()
			self.currentFile = cal
			self.currentFits = FakeFITS(self.currentFile)
			return True
		return False
	def getImage(self,extn):
		if self.currentFits is None:
			if self.allowMissing:
				return None
			else:
				raise ValueError
		return self.currentFits[extn]
	def getFileName(self):
		return self.currentFile

class FringeMap(CalibratorMap):
	'''Special case of CalibratorMap -- fringe images need be scaled to
	   match input images, but useful to keep the scaling mask locally'''
	def __init__(self,obsDb,calTab,nameMap=None,maskMap=None,sigThresh=1.0):
		super(FringeMap,self).__init__(obsDb,calTab,nameMap=nameMap,
		                               allowMissing=True)
		self.fringeMask = {}
		self.statsReg = stats_region(None)#'ccd_central_quadrant')
		self.sigThresh = sigThresh
		self.maskMap = maskMap
	def setTarget(self,f):
		changed = super(FringeMap,self).setTarget(f)
		if changed:
			if self.maskMap:
				self.maskFits = FakeFITS(self.maskMap(f))
			for extn in ['CCD%d' % i for i in range(1,5)]:
				im = self.currentFits[extn]
				mn,sig = array_stats(im[self.statsReg],method='median',
				                     rms=True,clip=True,
				                     clip_sig=5.0,clip_iters=1)
				self.fringeMask[extn] = np.abs((im-mn)/sig) < self.sigThresh
				self.fringeMask[extn] |= load_mask(self.maskFits[extn],
				                                   'nonzero')
	def getFringeScale(self,extn,inputIm):
		fringeIm = np.ma.array(self.currentFits[extn],
		                       mask=self.fringeMask[extn])
		# yet another sky calculation!
		skyVal = array_stats(inputIm[self.statsReg],method='median',clip=True)
		scaleIm = (inputIm-skyVal) / fringeIm
		scaleVal = array_stats(scaleIm,method='median',clip=True)
		return scaleVal

##############################################################################
#                                                                            #
# BokDataManager                                                             #
#   starts with an observations log and determines where to keep files, how  #
#   to map filenames, and provides an iterator over slices of the data.      #
#                                                                            #
##############################################################################

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
		self.fileFilter = None
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
		self.calDbFile = os.path.join(self.calDir,'caldb.pkl')
		self.calMap = {}
		try:
			self.calDb = load_caldb(self.calDbFile)
			self._config_cals()
		except:
			self.calDb = None
		# default filters for fringe corr
		self.fringeFilt = ['r','bokr','i','z'] 
		if 'good' not in self.obsDb.colnames:
			self.obsDb['good'] = np.ones(len(self.obsDb),dtype=bool)
	def _config_cals(self):
		self.calNameMap = SimpleFileNameMap(None,self.calDir)
		# should this be configurable by the user?
		self.calMaskMaps = {'fringe':self('imgmask')}
		self.calTable = {}
		self.calTable['bias'] = Table(rows=[d[:3] for d in self.calDb['zero']],
		                              names=('fileName','utDate','mjd'))
		cols = ('fileName','utDate','mjd','filter')
		for calType in ['bias','flat','illum','fringe','skyflat']:
			n,k = (3,'zero') if calType == 'bias' else (4,calType)
			try:
				dat = self.calDb[k]
				self.calTable[calType] = Table(rows=[d[:n] for d in dat],
				                               names=cols[:n])
				self.setCalMap(calType,'mjd',
				               maskMap=self.calMaskMaps.get(calType))
			except KeyError:
				pass
	def initCalDb(self):
		if not self.calDb:
			self.calDb = init_cal_db(self.calDbFile,self.obsDb,self.allFilt)
			self._config_cals()
	def setProcessSteps(self,steps):
		# XXX annoyingly hacky way to track where we're at in processing...
		self._all_process_steps = steps
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
	def getTmpDir(self):
		return self._tmpDir
	def setUtDates(self,utDates):
		# allows a range of dates to be specified by simply shortening the
		# string, i.e., '2014' is equivalent to '2014*'
		# use set to make sure the list is unique
		utdlist = list(set([ utd for f in utDates 
		                           for utd in self.allUtDates 
		                             if utd.startswith(f) ]))
		self.utDates = sorted(utdlist)
	def getUtDates(self):
		return self.utDates
	def iterUtDates(self):
		for u in self.utDates:
			self._curUtDate = u
			yield u
		self._curUtDate = None
	def setUtDirs(self,utDirs):
		allUtDirs = self.getUtDirs()
		# see setUtDates()
		dirlist = list(set([ dir for f in utDirs 
		                           for dir in allUtDirs 
		                             if dir.startswith(f) ]))
		# translate the nightly directory list back to utdates
		utdlist = [ self.obsDb['utDate'][self.obsDb['utDir']==d] 
		              for d in dirlist ]
		utdlist = np.unique(np.concatenate(utdlist))
		self.utDates = sorted(utdlist)
	def getUtDirs(self):
		return sorted(np.unique(self.obsDb['utDir']))
	def setFilters(self,filt):
		self.filt = filt
	def getFilters(self):
		return self.filt
	def setFringeFilters(self,filt):
		self.fringeFilt = filt
	def getFringeFilters(self):
		# recompose the list in case self.filt is reduced from full list
		return [b for b in self.filt if b in self.fringeFilt]
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
	def setFileFilter(self,fileFilter):
		self.fileFilter = fileFilter
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
	def getCalSequences(self,calType):
		return [ (fn,
		          [os.path.join(self.obsDb['utDir'][i],
		                        self.obsDb['fileName'][i]) 
		              for i in seq if self.obsDb['good'][i]])
		            for fn,utd,mjd,filt,seq in self.calDb[calType] 
		              if utd in self.utDates and 
		                 ((calType=='zero') or (filt in self.filt)) ]
	def updateCalSequences(self,calType,status):
		failed = [ fn for fn,success in status if not success ]
		newList = list(filter(lambda tup: tup[0] not in failed, 
		                      self.calDb[calType]))
		self.calDb[calType] = newList
		write_caldb(self.calDbFile,self.calDb)
		self._config_cals()
	def storeCalibrator(self,calType,frames,useFilt=True):
		calFn = caldb_store(self.calDbFile,self.obsDb,calType,[frames],
		                    useFilt=useFilt)[0]
		# have to reload the database
		self.calDb = load_caldb(self.calDbFile)
		# and add the mapping for this calib
		self.calTable[calType] = Table(rows=[d[:4] 
		                                       for d in self.calDb[calType]],
		                               names=('fileName','utDate',
		                                      'mjd','filter'))
		self.setCalMap(calType,'mjd')
		return self.calNameMap(calFn)
	def setCalMap(self,calType,mapType,fileName=None,maskMap=None):
		if mapType == None:
			self.calMap[calType] = NullCalibrator()
		elif mapType == 'master':
			self.calMap[calType] = MasterCalibrator(self.calNameMap(fileName))
		elif mapType == 'mjd':
			if calType == 'fringe':
				self.calMap[calType] = FringeMap(self.obsDb,
				                                 self.calTable[calType],
				                                 self.calNameMap,
				                                 maskMap=maskMap)
			else:
				self.calMap[calType] = CalibratorMap(self.obsDb,
				                                     self.calTable[calType],
				                                     self.calNameMap)
		else:
			raise ValueError
	def getCalMap(self,calType):
		return self.calMap[calType]
	def getCalPath(self,calFn):
		return os.path.join(self.calDir,calFn)
	def setInPlace(self,inPlace):
		self.inPlace = inPlace
		if self.inPlace:
			self.filesToMap = ['oscan','pass1cat','weight','imgmask',
			                   'skymask','skyfit','wcscat','cat','psf']
		else:
			self.filesToMap = self.fileSuffixes.keys()
	def __call__(self,t):
		procDir = self.procDir
		if t in self._all_process_steps:
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
		if t == 'raw':
			# special case to handle moving files from raw directory
			return SimpleFileNameMap(self.rawDir,procDir,fromRaw=True)
		elif t == 'cal':
			return self.calNameMap
		elif t in self.filesToMap:
			return SimpleFileNameMap(self.rawDir,procDir,
			                         self.fileSuffixes[t])
		else:
			return SimpleFileNameMap(self.rawDir,procDir)
	def getFiles(self,imType=None,utd=None,filt=None,
	             im_range=None,filterFun=None,includebad=False,
	             with_objnames=False,with_frames=False,as_sequences=False):
		if includebad:
			file_sel = np.ones_like(self.obsDb['good'])
		else:
			file_sel = self.obsDb['good'].copy()
		# select on UT date(s)
		if utd is not None:
			if isinstance(utd,basestring):
				utds = [utd]
			else:
				utds = utd
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
		if file_sel.any() and (filterFun or self.fileFilter):
			if filterFun is None:
				filterFun = self.fileFilter
			ii = np.where(file_sel)[0]
			file_sel[ii] &= filterFun(self.obsDb,ii)
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

