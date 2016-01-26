#!/usr/bin/env python

import os
import shutil
import numpy as np
import fitsio

from . import bokutil
from .bokproc import BokImArith

def smooth_correction_im(gradientFile,outputFile):
	from scipy.interpolate import LSQBivariateSpline
	spOrder = {'IM1':1,'IM2':1,'IM3':1,'IM4':1,
	           'IM5':3,'IM6':3,'IM7':1,'IM8':1,
	           'IM9':3,'IM10':3,'IM11':3,'IM12':3,
	           'IM13':1,'IM14':1,'IM15':1,'IM16':1}
	spKnots = {'IM1':2,'IM2':2,'IM3':2,'IM4':2,
	           'IM5':5,'IM6':5,'IM7':2,'IM8':2,
	           'IM9':7,'IM10':7,'IM11':5,'IM12':5,
	           'IM13':2,'IM14':2,'IM15':2,'IM16':2}
	gradientFits = fitsio.FITS(gradientFile)
	if os.path.exists(outputFile):
		os.unlink(outputFile)
	corrFits = fitsio.FITS(outputFile,'rw')
	corrFits.write(None,header=gradientFits[0].read_header())
	for hdu in gradientFits[1:]:
		extn = hdu.get_extname().upper()
		print 'fitting ramp correction to ',extn
		gim = hdu.read()
		ny,nx = gim.shape
		yi,xi = np.indices(gim.shape)
		mask = gim == 0 # need to use BokProcess to get bad pixel mask
		ii = np.where(~mask)
		nKnots = spKnots[extn]
		order = spOrder[extn]
		xknots = np.linspace(0,nx,nKnots)
		yknots = np.linspace(0,ny,nKnots)
		spfit = LSQBivariateSpline(xi[ii],yi[ii],gim[ii],
		                           xknots,yknots,kx=order,ky=order)
		im = spfit(np.arange(nx),np.arange(ny)).T
		corrFits.write(im,extname=extn,header=hdu.read_header())
	gradientFits.close()
	corrFits.close()

def make_rampcorr_image(dataMap,**kwargs):
	imNum = 1
	tmpDir = os.path.join(dataMap._tmpDir,'biasramps')
	inputMap = dataMap('oscan')
	if not os.path.exists(tmpDir):
		os.makedirs(tmpDir)
	rampFiles = []
	tmpRampFile = os.path.join(dataMap._tmpDir,'tmpbiasramp.fits')
	rampFile = os.path.join(dataMap.calDir,
	                        dataMap.getMaster('BiasRamp',name=True))
	for utd in dataMap.iterUtDates():
		files,frames = dataMap.getFiles(imType='zero',with_frames=True)
		if files is None:
			continue
		splits = np.where(np.diff(frames)>1)[0]
		if len(splits)==0:
			bias_seqs = [ files ]
		else:
			bias_seqs = np.split(files,splits+1)
		for biasNum,biasFiles in enumerate(bias_seqs,start=1):
			if len(biasFiles) < 2:
				continue
			outFile = os.path.join(tmpDir,'bias_%d_%s.fits' % (imNum,utd))
			imNum += 1
			rampFiles.append(outFile)
			if os.path.exists(outFile):
				continue
			imsub = BokImArith('-',inputMap(files[1]),
			                   output_map=lambda f: outFile)
			imsub.process_files([inputMap(files[0])])
	stackFun = bokutil.ClippedMeanStack()
	stackFun.stack(rampFiles,tmpRampFile)
	if kwargs.get('dosplinefit',False):
		smooth_correction_im(tmpRampFile,rampFile)
	else:
		shutil.copy(tmpRampFile,rampFile)

