#!/usr/bin/env python

import os
import re
import glob
import shutil
from copy import copy
from functools import partial
import multiprocessing
import numpy as np
from scipy.ndimage.morphology import binary_dilation
import fitsio
from astropy.table import Table

from .bokoscan import BokOverscanSubtract,BokOverscanSubtractWithSatFix
from . import bokio
from . import bokutil
from . import bokproc
from .bokdm import SimpleFileNameMap,BokDataManager
from . import bokastrom
from . import bokphot
from . import bokgnostic
from . import bokmkimage

all_process_steps = ['oscan','bias2d','flat2d','ramp',
                     'proc1','comb','fringe','skyflat',
                     'proc2','sky','wcs','cat']

def makeccd4image(dataMap,inputFile,outputFile=None,**kwargs):
	try:
		kwargs.pop('procmap') # only doing single image anyway
	except:
		pass
	if outputFile is None:
		ccd4map = bokio.FileNameMap(dataMap.getCalDir(),'_4ccd')
	else:
		ccd4map = lambda f: outputFile
	bokproc.combine_ccds([inputFile,],output_map=ccd4map,**kwargs)

def overscan_subtract(dataMap,fixsaturation=False,**kwargs):
	if fixsaturation:
		oscanSubtract = BokOverscanSubtractWithSatFix(input_map=dataMap('raw'),
	                                        output_map=dataMap('oscan'),
	                                        **kwargs)
	else:
		oscanSubtract = BokOverscanSubtract(input_map=dataMap('raw'),
	                                        output_map=dataMap('oscan'),
	                                        **kwargs)
	oscanSubtract.process_files(dataMap.getFiles())

def _bias_worker(dataMap,biasStack,nSkip,writeccdim,biasIn,**kwargs):
	biasFile,biasList = biasIn
	biasStack.stack(biasList[nSkip:],biasFile)
	if kwargs.get('verbose',0) >= 1:
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] 2DBIAS: %s' % (pid,biasFile)
		print '\n'.join(biasList[nSkip:])
	if writeccdim:
		makeccd4image(dataMap,biasFile,**kwargs)

def _bias_worker_exc(*args,**kwargs):
	try:
		_bias_worker(*args,**kwargs)
	except Exception,e:
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] 2DBIAS: %s FAILED!!! [%s]' % (pid,args[-1][0],e)
		return (args[-1][0],False)
	return (args[-1][0],True)

def make_2d_biases(dataMap,nSkip=2,reject='sigma_clip',
                   writeccdim=False,**kwargs):
	# need to make sure num processes is reset to 1 before calling these
	# routines since they will execute from within a subprocess
	procmap = kwargs.pop('procmap')
	_kwargs = copy(kwargs)
	_kwargs['processes'] = 1
	_kwargs['procmap'] = map
	biasStack = bokproc.BokBiasStack(input_map=dataMap('oscan'),
	                                 output_map=dataMap('cal'),
	                                 reject=reject,
	                                 find_rolloffs=True,
                                     **_kwargs)
	p_bias_worker = partial(_bias_worker_exc,dataMap,biasStack,
	                        nSkip,writeccdim,**kwargs)
	# returns [(biasFile,biasFiles)...]
	status = procmap(p_bias_worker,dataMap.getCalSequences('zero'))
	dataMap.updateCalSequences('zero',status)

def _flat_worker(dataMap,bias2Dsub,flatStack,normFlat,nSkip,writeccdim,
                 debug,flatIn,**kwargs):
	flatFile,flatList = flatIn
	if bias2Dsub:
		bias2Dsub.process_files(flatList)
	flatStack.stack(flatList[nSkip:],flatFile)
	if normFlat:
		raise NotImplementedError # currently all name mappings are broken
		if debug:
			shutil.copy(flatFile,
			            flatFile.replace('.fits','_raw.fits'))
		normFlat.process_files([flatFile])
	if kwargs.get('verbose',0) >= 1:
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] DOMEFLAT: %s' % (pid,flatFile)
		print '\n'.join(flatList[nSkip:])
	if writeccdim:
		makeccd4image(dataMap,flatFile,**kwargs)

def _flat_worker_exc(*args,**kwargs):
	flatFile = args[-1][0]
	try:
		_flat_worker(*args,**kwargs)
	except Exception,e:
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] DOMEFLAT: %s FAILED!!! [%s]' % (pid,flatFile,e)
		flatStack = args[2]
		flatFilePath = flatStack.outputNameMap(flatFile)
		try:
			os.remove(flatFilePath)
		except:
			pass
		return (flatFile,False)
	return (flatFile,True)

def make_dome_flats(dataMap,nobiascorr=False,
                    nSkip=1,reject='sigma_clip',writeccdim=False,
	                usepixflat=True,debug=False,**kwargs):
	# need to make sure num processes is reset to 1 before calling these
	# routines since they will execute from within a subprocess
	procmap = kwargs.pop('procmap')
	_kwargs = copy(kwargs)
	_kwargs['processes'] = 1
	_kwargs['procmap'] = map
	if nobiascorr:
		bias2Dsub = None
	else:
		bias2Dsub = bokproc.BokCCDProcess(input_map=dataMap('oscan'),
		                                  output_map=dataMap('bias'),
		                                  bias=dataMap.getCalMap('bias'),
		                                  header_key='BIAS2D',
		                                  **_kwargs)
	flatStack = bokproc.BokDomeFlatStack(reject=reject,
	                                     input_map=dataMap('bias'),
	                                     output_map=dataMap('cal'),
	                                     **_kwargs)
	if usepixflat:
		if debug:
			ffmap = SimpleFileNameMap(None,dataMap.procDir,'_fit')
			bfmap = SimpleFileNameMap(None,dataMap.procDir,'_binned')
		else:
			ffmap,bfmap = None,None
		normFlat = bokproc.NormalizeFlat(_normed_flat_fit_map=ffmap,
		                                 _binned_flat_map=bfmap,**kwargs)
	else:
		normFlat = None
	p_flat_worker = partial(_flat_worker_exc,dataMap,bias2Dsub,flatStack,
	                        normFlat,nSkip,writeccdim,debug,**kwargs)
	status = procmap(p_flat_worker,dataMap.getCalSequences('flat'))
	dataMap.updateCalSequences('flat',status)

def _ramp_worker(tmpDir,inputMap,calMap,verbose,biasIn):
	biasFile,biasList = biasIn
	if verbose >= 1:
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] RAMP: %s' % (pid,biasFile)
	outFile = os.path.join(tmpDir,'ramp'+biasFile+'.fits')
	if os.path.exists(outFile):
		return outFile
	imsub = bokproc.BokImArith('-',calMap(biasFile),
	                           output_map=lambda f: outFile)
	imsub.process_files([inputMap(biasList[0])])
	return outFile

# XXX removed smoothing of image b/c not being used, but could do with spline smoother?
def make_rampcorr_image(dataMap,**kwargs):
	processes = kwargs.get('processes',1)
	procmap = kwargs.pop('procmap')
	tmpDir = os.path.join(dataMap._tmpDir,'biasramps')
	inputMap = dataMap('oscan')
	calMap = dataMap('cal')
	if not os.path.exists(tmpDir):
		os.makedirs(tmpDir)
	rampFile = dataMap.getCalMap('ramp').getFileName()
	p_ramp_worker = partial(_ramp_worker,tmpDir,inputMap,calMap,
	                        kwargs.get('verbose',0))
	rampFiles = procmap(p_ramp_worker,dataMap.getCalSequences('zero'))
	stackFun = bokutil.ClippedMeanStack()
	stackFun.stack(rampFiles,rampFile)

def balance_gains(dataMap,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=dataMap('proc1'),
	                                     mask_map=dataMap.getCalMap('badpix'),
	                                                **kwargs)
	gainMap = {'corrections':{},'skyvals':{}}
	for utd in dataMap.iterUtDates():
		files,ii = dataMap.getFiles(imType='object',with_frames=True)
		if files is None:
			continue
		filt = dataMap.obsDb['filter'][ii]
		diagfile = os.path.join(dataMap.getDiagDir(), 'gainbal_%s.npz'%utd)
		if os.path.exists(diagfile):
			gainDat = np.load(diagfile)
			gainCor = gainDat['gainCor']
			skyV = gainDat['skys']
		else:
			gainBalance.process_files(files,filt)
			gainCor = gainBalance.calc_mean_corrections()
			ampGainV,ccdGainV,gainCorV,ampTrend,ccdTrend,skyV = \
			               gainBalance.get_values()
		for f,gc,skyv in zip(files,gainCor,skyV):
			gainMap['corrections'][f] = gc
			gainMap['skyvals'][f] = skyv
		gainBalance.reset()
		if not os.path.exists(diagfile) and \
		     not kwargs.get('nosavegain',False):
			np.savez(diagfile,files=[os.path.basename(f) for f in files],
			         gains=gainCorV,skys=skyV,gainCor=gainCor,
			         ampTrend=ampTrend,ccdTrend=ccdTrend,
			         rawAmpGain=ampGainV,rawCcdGain=ccdGainV)
	return gainMap

def files_by_utdfilt(dataMap,imType='object',filt=None):
	if len(dataMap.utDates) > 20: # XXX >> nProc
		if filt is None:
			filt = dataMap.getFilters()
		filesUtdFilt = [ dataMap.getFiles(imType=imType,filt=_filt)
		                   for _utd in dataMap.iterUtDates()
		                     for _filt in filt ]
		filesUtdFilt = filter(lambda l: l is not None, filesUtdFilt)
		# also return flattened list
		files = [ f for utdfilt in filesUtdFilt for f in utdfilt ]
	else:
		files = filesUtdFilt = dataMap.getFiles(imType=imType,filt=filt)
	return files,filesUtdFilt

def process_all(dataMap,nobiascorr=False,noflatcorr=False,
                fixpix=False,rampcorr=False,noweightmap=False,
                nocombine=False,prockey='CCDPROC',**kwargs):
	# 1. basic processing (bias and flat-field correction, fixpix, 
	#    nominal gain correction
	bias = None if nobiascorr else dataMap.getCalMap('bias')
	flat = None if noflatcorr else dataMap.getCalMap('flat')
	ramp = None if not rampcorr else dataMap.getCalMap('ramp')
	proc = bokproc.BokCCDProcess(input_map=dataMap('oscan'),
	                             output_map=dataMap('proc1'),
	                             mask_map=dataMap.getCalMap('badpix'),
	                             header_key=prockey,
	                             bias=bias,flat=flat,ramp=ramp,
	                             fringe=None,illum=None,skyflat=None,
	                             fixpix=fixpix,**kwargs)
	files,filesUtdFilt = files_by_utdfilt(dataMap)
	if files is None or len(files)==0:
		return
	proc.process_files(filesUtdFilt)
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
	if not noweightmap:
		# 4. construct weight maps starting from raw images
		whmap = bokproc.BokWeightMap(input_map=dataMap('raw'),
		                             output_map=dataMap('weight'),
		                             flat=flat,
		                             _mask_map=dataMap.getCalMap('badpix'),
		                             **kwargs)
		whmap.process_files(filesUtdFilt)
		# rescale the gain corrections to inverse variance
		for f in gainMap['corrections']:
			gainMap['corrections'][f] = (gainMap['corrections'][f].copy())**-2
		bokproc.combine_ccds(files,
		                     input_map=dataMap('weight'), 
		                     output_map=dataMap('weight'), 
		                     gain_map=gainMap,
		                     **kwargs)

def make_illumcorr_image(dataMap,byUtd=True,filterFun=None,
                         min_images=10,max_images=None,
                         iterfit=True,**kwargs):
	redo = kwargs.get('redo',False)
	verbose = kwargs.get('verbose',0)
	if byUtd:
		filtAndUtd = [ (f,u) for f in dataMap.getFilters()
		                       for u in dataMap.getUtDates() ]
	else:
		filtAndUtd = [ (f,None) for f in dataMap.getFilters()]
	for filt,utd in filtAndUtd:
		files,frames = dataMap.getFiles(imType='object',filt=filt,utd=utd,
		                                filterFun=filterFun,
		                                with_frames=True)
		if files is None or len(files) < min_images:
			continue
		outFn = dataMap.storeCalibrator('illum',frames)
		if os.path.exists(outFn) and not redo:
			if verbose > 0:
				print 'ILLUM: %s already exists' % outFn
			continue
		if verbose > 0:
			print 'ILLUM: generating %s from %d images' % (outFn,len(files))
		#
		tmpFn = 'tmp'+os.path.basename(outFn)
		tmpSkyFlatFile = os.path.join(dataMap._tmpDir,tmpFn)
		stackFun = bokutil.ClippedMeanStack(input_map=dataMap('comb'),
		                                    scale='normalize_median',
		                                    clip_iters=3,clip_sig=2.0,
		                                    **kwargs)
		if max_images is not None:
			# keep the image list in the same order
			ii = np.random.randint(0,max_images)
			files = [files[i] for i in ii]
		print 'stacking %d files for illumination' % (len(files))
		stackFun.stack(files,tmpSkyFlatFile)
		fits = bokutil.BokMefImage(tmpSkyFlatFile,
		                           mask_file=dataMap.getCalMap('badpix4'),
		                           read_only=True)
		illum = bokproc.SplineBackgroundFit(fits,nKnots=7,order=3,nbin=8)
		if iterfit:
			mask = {}
			for extn,data,hdr in fits:
				resid = data - illum.get(extn)
				arr = bokutil.array_clip(resid[10:-10:10,10:-10:10])
				mn = arr.mean()
				sd = arr.std()
				mask[extn] = np.abs(resid-mn) > 3*sd
			fits.add_mask(mask)
			illum = bokproc.SplineBackgroundFit(fits,nKnots=17,order=3,nbin=4)
		normFun = lambda arr: arr / np.float32(illum(0,0))
		illum.write(outFn,opfun=normFun,clobber=True)

def make_fringe_masters(dataMap,byUtd=False,**kwargs):
	caldir = dataMap.getCalDir()
	stackin = dataMap('fringe') # XXX
	fringeStack = bokproc.BokFringePatternStack(input_map=stackin,
#	                                            mask_map=dataMap('skymask'),
	                       raw_stack_file=bokio.FileNameMap(caldir,'_raw'),
	                                        header_bad_key='BADSKY',
	                                            **kwargs)
	fringeStack.set_badpixelmask(dataMap.getCalMap('badpix4'))
	if byUtd:
		filtAndUtd = [ (f,u) for f in dataMap.getFringeFilters()
		                       for u in dataMap.getUtDates() ]
	else:
		filtAndUtd = [ (f,None) for f in dataMap.getFringeFilters()]
	for filt,utd in filtAndUtd:
		files,frames = dataMap.getFiles(imType='object',filt=filt,utd=utd,
		                                with_frames=True)
		outfn = dataMap.storeCalibrator('fringe',frames)
		fringeStack.stack(files,outfn)

def make_supersky_flats(dataMap,byUtd=False,interpFill=True,**kwargs):
	caldir = dataMap.getCalDir()
	stackin = dataMap('sky') # XXX
	statsReg = bokutil.stats_region(None,16)
	growKern = None
	growThresh = 2.0
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=dataMap('skymask'),
	                    exposure_time_map=bokio.FileNameMap(caldir,'.exp'),
	                                        header_bad_key='BADSKY',
	                                            **kwargs)
	skyFlatStack.set_badpixelmask(dataMap.getCalMap('badpix4'))
	if byUtd:
		filtAndUtd = [ (f,u) for f in dataMap.getFilters()
		                       for u in dataMap.getUtDates() ]
	else:
		filtAndUtd = [ (f,None) for f in dataMap.getFilters()]
	for filt,utd in filtAndUtd:
		files,frames = dataMap.getFiles(imType='object',filt=filt,utd=utd,
		                                with_frames=True)
		_outfn = outfn = dataMap.storeCalibrator('skyflat',frames)
		if interpFill:
			_outfn = _outfn.replace('.fits','_raw.fits')
		skyFlatStack.stack(files,_outfn)
		if interpFill:
			fits = bokutil.BokMefImage(_outfn,
			                           output_file=outfn,
			                           mask_file=dataMap.getCalMap('badpix4'))
			expMapFn = bokio.FileNameMap(caldir,'.exp')(_outfn)
			expFits = bokutil.FakeFITS(expMapFn)
			mask = {}
			for ccd,data,hdr in fits:
				expim = expFits[ccd]
				# mask underexposed regions
				expmask = expim < 0.5*np.median(expim[expim>0])
				# mask large deviations
				statsPix = bokutil.array_clip(data[statsReg],
				                              clip_iters=2,clip_sig=5.0)
				snr = np.abs(data-statsPix.mean())/statsPix.std()
				m = (snr > 10) | (data==1) | expmask
				# grow mask
				binary_dilation(m,mask=(snr>growThresh),iterations=0,
				                structure=growKern,output=m)
				mask[ccd] = m
			fits.add_mask(mask)
			backfit = bokproc.SplineBackgroundFit(fits,nKnots=50,
			                                      order=1,nbin=16)
			for ccd,data,hdr in fits:
				data[data.mask] = backfit.get(ccd)[data.mask]
				fits.update(data.filled(),hdr)
			fits.close()

def process_all2(dataMap,skyArgs,noillumcorr=False,noskyflatcorr=False,
                 nofringecorr=False,noskysub=False,noweightmap=False,
                 prockey='CCDPRO2',redoskymask=False,save_sky=False,**kwargs):
	#
	# Second round: illumination, fringe, and skyflat corrections
	#
	fringe = None if nofringecorr else dataMap.getCalMap('fringe')
	illum = None if noillumcorr else dataMap.getCalMap('illum')
	skyflat = None if noskyflatcorr else dataMap.getCalMap('skyflat')
	# the full list of files to process
	files,filesUtdFilt = files_by_utdfilt(dataMap)
	if files is None or len(files)==0:
		return
	# if only applying a fringe correction, only CCDProcess the images 
	# that need it
	if noillumcorr and noskyflatcorr and not nofringecorr:
		_,_filesUtdFilt = files_by_utdfilt(dataMap,
		                                 filt=dataMap.getFringeFilters())
	else:
		_filesUtdFilt = filesUtdFilt
	if len(_filesUtdFilt) > 0:
		proc = bokproc.BokCCDProcess(input_map=dataMap('comb'),
		                             output_map=dataMap('proc2'),
		                             mask_map=dataMap.getCalMap('badpix4'),
		                             header_key=prockey,
		                             gain_multiply=False,bias=None,flat=None,
		                             ramp=None,fixpix=False,fringe=fringe,
		                             illum=illum,skyflat=skyflat,
		                             divide_exptime=True,
		                             **kwargs)
		proc.process_files(_filesUtdFilt)
	if not noweightmap and not (noillumcorr and noskyflatcorr):
		# need to process the weight maps in the same fashion
		wtproc = bokproc.BokCCDProcess(input_map=dataMap('weight'), 
		                               output_map=dataMap('weight'), 
		                               mask_map=dataMap.getCalMap('badpix4'),
		                               header_key=prockey,
		                               gain_multiply=False,bias=None,flat=None,
		                               ramp=None,fixpix=False,fringe=None,
		                               illum=illum,skyflat=skyflat,
		                               divide_exptime=True,asweight=True,
		                               **kwargs)
		wtproc.process_files(filesUtdFilt)
	if noskysub:
		return
	#
	# Sky subtraction
	#
	# Generate sky masks by agressively masking objects
	skymskkwargs = copy(kwargs) # sky masks are time-consuming so only
	skymskkwargs['clobber'] = redoskymask  # redo if really necessary
	skyFlatMask = bokproc.BokGenerateSkyFlatMasks(
	                                    input_map=dataMap('proc2'),
	                                    output_map=dataMap('skymask'),
	                                    mask_map=dataMap.getCalMap('badpix4'),
	                                              **skymskkwargs)
	skyFlatMask.process_files(files)
	skyfitmap = dataMap('skyfit') if save_sky else None
	skySub = bokproc.BokSkySubtract(input_map=dataMap('proc2'),
	                                output_map=dataMap('sky'),
	                                mask_map=dataMap('skymask'),
	                                skyfit_map=skyfitmap,
	                                **dict(skyArgs.items()+kwargs.items()))
	skySub.add_mask(dataMap.getCalMap('badpix4'))
	skySub.process_files(files)

def _wcs_worker(dataMap,inputType,redowcscat,savewcs,keepwcscat,
                clobber,verbose,inp):
	try:
		imFile,fieldName = inp
		imageFile = dataMap(inputType)(imFile)
		catFile = dataMap('wcscat')(imFile)
		try:
			pid = multiprocessing.current_process().name.split('-')[1]
		except:
			pid = '1'
		print '[%2s] WCS: %s' % (pid,imFile)
		# kwargs sent to the following are added sextractor/scamp parameters
		#  (e.g., {'VERBOSE':'FULL'}), so remap the pipeline kwargs 
		bokphot.sextract(imageFile,catFile,
		                 clobber=redowcscat,verbose=verbose)
		bokastrom.scamp_solve(imageFile,catFile,
		                      dataMap.getScampRefCat(fieldName),
		                      filt='r',savewcs=savewcs,
		                      clobber=clobber,verbose=verbose)
		if not keepwcscat:
			os.unlink(catFile)
	except:
		pass

def set_wcs(dataMap,inputType='sky',savewcs=False,keepwcscat=True,
            redowcscat=False,**kwargs):
	procmap = kwargs.pop('procmap',map)
	filesAndFields = dataMap.getFiles(imType='object',with_objnames=True)
	p_wcs_worker = partial(_wcs_worker,dataMap,inputType,redowcscat,
	                       savewcs,keepwcscat,kwargs.get('clobber',False),
	                       kwargs.get('verbose',0))
	status = procmap(p_wcs_worker,zip(*filesAndFields))

def _cat_worker(dataMap,inputType,imFile,**kwargs):
	try:
		pid = multiprocessing.current_process().name.split('-')[1]
	except:
		pid = '1'
	print '[%2s] CATALOG: %s' % (pid,imFile)
	try:
		imageFile = dataMap(inputType)(imFile)
		psfFile = dataMap('psf')(imFile)
		if not os.path.exists(psfFile):
			catFile = dataMap('wcscat')(imFile)
			bokphot.sextract(imageFile,catFile,full=False,**kwargs)
			bokphot.run_psfex(catFile,psfFile,**kwargs)
		catFile = dataMap('cat')(imFile)
		bokphot.sextract(imageFile,catFile,psfFile,full=True,**kwargs)
	except:
		pass

def make_catalogs(dataMap,inputType='sky',**kwargs):
	procmap = kwargs.pop('procmap',map)
	files = dataMap.getFiles(imType='object')
	p_cat_worker = partial(_cat_worker,dataMap,inputType,**kwargs)
	status = procmap(p_cat_worker,files)

def bokpipe(dataMap,**kwargs):
	redo = kwargs.get('redo',False)
	steps = kwargs.get('steps')
	debug = kwargs.get('debug',False)
	verbose = kwargs.get('verbose',0)
	processes = kwargs.get('processes',1)
	procmap = kwargs.get('procmap')
	maxmem = kwargs.get('maxmem',5)
	chunkSize = 10
	if processes > 1:
		pool = multiprocessing.Pool(processes)
		# need to test this more to see if it improves efficiency, but
		# probably need to wrap it in an object
		#procmap = partial(pool.map,chunksize=chunkSize)
		procmap = pool.map
	else:
		procmap = map
	pipekwargs = {'clobber':redo,'verbose':verbose,
	              'processes':processes,'procmap':procmap,'maxmem':maxmem}
	# fixpix is sticking nan's into the images in unmasked pixels (???)
	fixpix = False #True
	writeccdims = kwargs.get('calccdims',False)
	timerLog = bokutil.TimerLog()
	biasMap = None
	if 'oscan' in steps:
		overscan_subtract(dataMap,fixsaturation=kwargs.get('fixsaturation'),
		                  **pipekwargs)
		timerLog('overscans')
	if 'bias2d' in steps:
		make_2d_biases(dataMap,writeccdim=writeccdims,**pipekwargs)
		timerLog('2d biases')
	if 'flat2d' in steps:
		make_dome_flats(dataMap,writeccdim=writeccdims,
		                nobiascorr=kwargs.get('nobiascorr',False),
		                usepixflat=not kwargs.get('nousepixflat',False),
		                debug=debug,**pipekwargs)
		timerLog('dome flats')
	if 'ramp' in steps:
		make_rampcorr_image(dataMap,**pipekwargs)
		timerLog('ramp correction')
	if 'proc1' in steps or 'comb' in steps:
		process_all(dataMap,
		            nobiascorr=kwargs.get('nobiascorr',False),
		            noflatcorr=kwargs.get('noflatcorr',False),
		            fixpix=fixpix,
		            rampcorr=kwargs.get('rampcorr',False),
		            nocombine=kwargs.get('nocombine',False),
		            gain_multiply=not kwargs.get('nogainmul',False),
		            nosavegain=kwargs.get('nosavegain',False),
		            noweightmap=kwargs.get('noweightmap',False),
		            prockey=kwargs.get('prockey','CCDPROC'),
		            **pipekwargs)
		timerLog('ccdproc')
	if 'illum' in steps:
		make_illumcorr_image(dataMap,
		                     filterFun=kwargs.get('illum_filter_fun'),
		                     **pipekwargs)
		timerLog('illumination corr')
	if 'fringe' in steps:
		make_fringe_masters(dataMap,**pipekwargs)
		timerLog('fringe masters')
	if 'skyflat' in steps:
		make_supersky_flats(dataMap,**pipekwargs)
		timerLog('supersky flats')
	if 'proc2' in steps:
		skyArgs = { k.lstrip('sky'):kwargs[k] 
		                 for k in ['skymethod','skyorder']}
		process_all2(dataMap,skyArgs,
		             noillumcorr=kwargs.get('noillumcorr'),
		             nofringecorr=kwargs.get('nofringecorr'),
		             noskyflatcorr=kwargs.get('noskyflatcorr'),
		             noweightmap=kwargs.get('noweightmap'),
		             noskysub=kwargs.get('noskysub'),
		             prockey=kwargs.get('prockey','CCDPRO2'),
		             redoskymask=kwargs.get('redoskymask'),
		             save_sky=kwargs.get('savesky'),
		             **pipekwargs)
		timerLog('process2')
	if 'wcs' in steps:
		set_wcs(dataMap,
		        savewcs=kwargs.get('savewcs',False),
		        keepwcscat=kwargs.get('keepwcscat',True),
		        redowcscat=kwargs.get('redowcscat',False),
		        **pipekwargs)
		timerLog('wcs')
	if 'cat' in steps:
		make_catalogs(dataMap,**pipekwargs)
		timerLog('catalog')
	timerLog.dump()
	if processes > 1:
		pool.close()

def _img_worker(imdir,_fmap,maskmap,fin):
	f = _fmap(fin)
	print f
	if not os.path.exists(f):
		return
	d,fn = os.path.split(f)
	d,utd = os.path.split(d)
	imgd = os.path.join(imdir,utd)
	imgfile = os.path.join(imgd,fn.replace('.fits','.png'))
	if os.path.exists(imgfile):
		return
	bokmkimage.make_fov_image_fromfile(f,imgfile,mask=maskmap(fin))

def make_images(dataMap,imtype='sky',msktype=None,processes=1):
	import matplotlib.pyplot as plt
	if processes > 1:
		pool = multiprocessing.Pool(processes)
		procmap = pool.map
	else:
		procmap = map
	_fmap = dataMap(imtype)
	if msktype=='badpix':
		msktype = 'badpix4' # XXX
	if msktype==None:
		maskmap = bokio.NullNameMap
	else:
		maskmap = dataMap.getCalMap(msktype)
	imdir = os.path.join(dataMap.getProcDir(),'images')
	for utDir in dataMap.getUtDirs():
		imgd = os.path.join(imdir,utDir)
		if not os.path.exists(imgd):
			os.makedirs(imgd)
	plt.ioff()
	p_img_worker = partial(_img_worker,imdir,_fmap,maskmap)
	procmap(p_img_worker,dataMap.getFiles(imType='object'))
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
	dataMap.setProcessSteps(all_process_steps) # XXX argh
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
		for d in [dataMap.procDir,dataMap.getCalDir(),dataMap.getDiagDir(),
		          dataMap.getTmpDir()]:
			if not os.path.exists(d):
				os.makedirs(d)
		for _utdir in dataMap.getUtDirs():
			utdir = os.path.join(dataMap.procDir,_utdir)
			if not os.path.exists(utdir): os.mkdir(utdir)
	#
	dataMap.initCalDb()
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
	parser.add_argument('--maxmem',type=float,default=5,
	                help='maximum memory in GB for stacking images')
	parser.add_argument('--calccdims',action='store_true',
	                help='generate CCD-combined images for calibration data')
	parser.add_argument('--fixsaturation',action='store_true',
	                help='correct overflowed pixels to have saturation value')
	parser.add_argument('--nobiascorr',action='store_true',
	                help='do not apply bias correction')
	parser.add_argument('--noflatcorr',action='store_true',
	                help='do not apply flat correction')
	parser.add_argument('--rampcorr',action='store_true',
	                help='apply bias ramp correction')
	parser.add_argument('--noillumcorr',action='store_true',
	                help='do not apply illumination correction')
	parser.add_argument('--nofringecorr',action='store_true',
	                help='do not apply fringe correction')
	parser.add_argument('--noskyflatcorr',action='store_true',
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
	parser.add_argument('--redoskymask',action='store_true',
	                help='redo sky mask generation')
	parser.add_argument('--savesky',action='store_true',
	                help='save sky background fit')
	parser.add_argument('--noweightmap',action='store_true',
	                help='do not generate weight maps')
	parser.add_argument('--prockey',type=str,default=None,
	                help='set new header key for ccdproc')
	parser.add_argument('--tmpdirin',action='store_true',
	                help='read files from temporary directory')
	parser.add_argument('--tmpdirout',action='store_true',
	                help='write files to temporary directory')
	parser.add_argument('--images',action='store_true',
	                help='make png images instead of processing ')
	parser.add_argument('--imagetype',type=str,default='sky',
	                help='make images from (imtype,[msktype]) [default: sky]')
	parser.add_argument('--savewcs',action='store_true',
	                help='write wcs to headers')
	parser.add_argument('--wcscheck',action='store_true',
	                help='make astrometry diagnostic files')
	return parser

def run_pipe(dataMap,args,**_kwargs):
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
	for k,v in _kwargs.items():
		kwargs[k] = v
	# run pipeline processes
	if args.images:
		make_images(dataMap,*args.imagetype.split(','),
		            processes=args.processes)
	elif args.wcscheck:
		files = [dataMap('sky')(f) for f in dataMap.getFiles('object')]
		bokgnostic.run_scamp_diag(files)
	else:
		bokpipe(dataMap,**kwargs)

