#!/usr/bin/env python

import os
import re
import glob
import shutil
from copy import copy
from functools import partial
import multiprocessing
import numpy as np
import fitsio
from astropy.table import Table

from .bokoscan import BokOverscanSubtract,BokOverscanSubtractWithSatFix
from . import bokio
from . import bokutil
from . import bokproc
from . import bokrampcorr
from . import bokillumcorr
from .bokdm import SimpleFileNameMap,BokDataManager
from . import bokastrom
from . import bokphot
from . import bokgnostic
from . import bokmkimage

all_process_steps = ['oscan','bias2d','flat2d',
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
                                     **_kwargs)
	p_bias_worker = partial(_bias_worker,dataMap,biasStack,
	                        nSkip,writeccdim,**kwargs)
	# returns [(biasFile,biasFiles)...]
	procmap(p_bias_worker,dataMap.getCalSequences('zero'))

def _flat_worker(dataMap,bias2Dsub,flatStack,normFlat,nSkip,writeccdim,
                 debug,flatIn,**kwargs):
	flatFile,flatList = flatIn
	if bias2Dsub:
		bias2Dsub.process_files(flatList)
	flatStack.stack(flatList[nSkip:],flatFile)
	if normFlat:
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
	p_flat_worker = partial(_flat_worker,dataMap,bias2Dsub,flatStack,
	                        normFlat,nSkip,writeccdim,debug,**kwargs)
	procmap(p_flat_worker,dataMap.getCalSequences('flat'))

def balance_gains(dataMap,**kwargs):
	# need bright star mask here?
	gainBalance = bokproc.BokCalcGainBalanceFactors(
	                                     input_map=dataMap('proc1'),
	                                     mask_map=dataMap.getCalMap('badpix'),
	                                                **kwargs)
	gainMap = {'corrections':{},'skyvals':{}}
	for utd in dataMap.iterUtDates():
		files = dataMap.getFiles(imType='object')
		if files is None:
			continue
		diagfile = os.path.join(dataMap.getDiagDir(), 'gainbal_%s.npz'%utd)
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
	if not noweightmap:
		# 4. construct weight maps starting from raw images
		whmap = bokproc.BokWeightMap(input_map=dataMap('raw'),
		                             output_map=dataMap('weight'),
		                             flat=flat,
		                             _mask_map=dataMap.getCalMap('badpix'),
		                             **kwargs)
		whmap.process_files(files)
		# rescale the gain corrections to inverse variance
		for f in gainMap['corrections']:
			gainMap['corrections'][f] = (gainMap['corrections'][f].copy())**-2
		bokproc.combine_ccds(files,
		                     input_map=dataMap('weight'), 
		                     output_map=dataMap('weight'), 
		                     gain_map=gainMap,
		                     **kwargs)

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
		filtAndUtd = [ fu for fu in zip(dataMap.getFringeFilters(),
		                                dataMap.getUtDates()) ]
	else:
		filtAndUtd = [ (f,None) for f in dataMap.getFringeFilters()]
	for filt,utd in filtAndUtd:
		files,frames = dataMap.getFiles(imType='object',filt=filt,utd=utd,
		                                with_frames=True)
		outfn = dataMap.storeCalibrator('fringe',frames)
		fringeStack.stack(files,outfn)

def make_supersky_flats(dataMap,byUtd=False,**kwargs):
	caldir = dataMap.getCalDir()
	stackin = dataMap('sky') # XXX
	skyFlatStack = bokproc.BokNightSkyFlatStack(input_map=stackin,
	                                            mask_map=dataMap('skymask'),
	                    exposure_time_map=bokio.FileNameMap(caldir,'.exp'),
	                       raw_stack_file=bokio.FileNameMap(caldir,'_raw'),
	                                        header_bad_key='BADSKY',
	                                            **kwargs)
	skyFlatStack.set_badpixelmask(dataMap.getCalMap('badpix4'))
	if byUtd:
		filtAndUtd = [ fu for fu in zip(dataMap.getFilters(),
		                                dataMap.getUtDates()) ]
	else:
		filtAndUtd = [ (f,None) for f in dataMap.getFilters()]
	for filt,utd in filtAndUtd:
		files,frames = dataMap.getFiles(imType='object',filt=filt,utd=utd,
		                                with_frames=True)
		outfn = dataMap.storeCalibrator('skyflat',frames)
		skyFlatStack.stack(files,outfn)

def process_all2(dataMap,skyArgs,noillumcorr=False,noskyflatcorr=False,
                 nofringecorr=False,noskysub=False,noweightmap=False,
                 prockey='CCDPRO2',save_sky=False,**kwargs):
	#
	# Second round: illumination, fringe, and skyflat corrections
	#
	fringe = None if nofringecorr else dataMap.getCalMap('fringe')
	illum = None if noillumcorr else dataMap.getCalMap('illum')
	skyflat = None if noskyflatcorr else dataMap.getCalMap('skyflat')
	# the full list of files to process
	files = dataMap.getFiles(imType='object')
	if files is None:
		return
	# if only applying a fringe correction, only CCDProcess the images 
	# that need it
	if noillumcorr and noskyflatcorr and not nofringecorr:
		_files = dataMap.getFiles(imType='object',
		                          filt=dataMap.getFringeFilters())
	else:
		_files = files
	if len(_files) > 0:
		proc = bokproc.BokCCDProcess(input_map=dataMap('comb'),
		                             output_map=dataMap('proc2'),
		                             mask_map=dataMap.getCalMap('badpix4'),
		                             header_key=prockey,
		                             gain_multiply=False,bias=None,flat=None,
		                             ramp=None,fixpix=False,fringe=fringe,
		                             illum=illum,skyflat=skyflat,
		                             **kwargs)
		proc.process_files(_files)
	if not noweightmap and not (noillumcorr and noskyflatcorr):
		# need to process the weight maps in the same fashion
		wtproc = bokproc.BokCCDProcess(input_map=dataMap('weight'), 
		                               output_map=dataMap('weight'), 
		                               mask_map=dataMap.getCalMap('badpix4'),
		                               header_key=prockey,
		                               gain_multiply=False,bias=None,flat=None,
		                               ramp=None,fixpix=False,fringe=None,
		                               illum=illum,skyflat=skyflat,
		                               asweight=True,
		                               **kwargs)
		wtproc.process_files(files)
	if noskysub:
		return
	#
	# Sky subtraction
	#
	# Generate sky masks by agressively masking objects
	skyFlatMask = bokproc.BokGenerateSkyFlatMasks(
	                                    input_map=dataMap('proc2'),
	                                    output_map=dataMap('skymask'),
	                                    mask_map=dataMap.getCalMap('badpix4'),
	                                              **kwargs)
	skyFlatMask.process_files(files)
	skyfitmap = dataMap('skyfit') if save_sky else None
	skySub = bokproc.BokSkySubtract(input_map=dataMap('proc2'),
	                                output_map=dataMap('sky'),
	                                mask_map=dataMap('skymask'),
	                                skyfit_map=skyfitmap,
	                                **dict(skyArgs.items()+kwargs.items()))
	skySub.add_mask(dataMap.getCalMap('badpix4'))
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

def bokpipe(dataMap,**kwargs):
	redo = kwargs.get('redo',False)
	steps = kwargs.get('steps')
	debug = kwargs.get('debug',False)
	verbose = kwargs.get('verbose',0)
	processes = kwargs.get('processes',1)
	procmap = kwargs.get('procmap')
	maxmem = kwargs.get('maxmem',5)
	if processes > 1:
		pool = multiprocessing.Pool(processes)
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
	if not os.path.exists(imgd):
		os.mkdir(imgd)
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
	if not os.path.exists(imdir):
		os.mkdir(imdir)
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
		for d in [dataMap.procDir,dataMap.getCalDir(),dataMap.getDiagDir()]:
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
	parser.add_argument('--savesky',action='store_true',
	                help='save sky background fit')
	parser.add_argument('--noweightmap',action='store_true',
	                help='do not generate weight maps')
	parser.add_argument('--prockey',type=str,default=None,
	                help='set new header key for ccdproc')
	parser.add_argument('--tmpdirin',action='store_true',
	                help='read files from temporary directory')
	parser.add_argument('--tmpdirout',
	                help='write files to temporary directory')
	parser.add_argument('--images',action='store_true',
	                help='make png images instead of processing ')
	parser.add_argument('--imagetype',type=str,default='sky',
	                help='make images from (imtype,[msktype]) [default: sky]')
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
		make_images(dataMap,*args.imagetype.split(','),
		            processes=args.processes)
	elif args.makerampcorr:
		bokrampcorr.make_rampcorr_image(dataMap)#,**kwargs)
	elif args.makeillumcorr:
		bokillumcorr.make_illumcorr_image(dataMap)#,**kwargs)
	elif args.wcscheck:
		files = [dataMap('sky')(f) for f in dataMap.getFiles('object')]
		bokgnostic.run_scamp_diag(files)
	else:
		bokpipe(dataMap,**kwargs)

