#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import LSQBivariateSpline
from astropy.modeling import models,fitting

from bokpipe import *

def make_skyflat(skyFlatFile):
	'''Make a dark sky flat template for the illumination correction'''
	stackPars = {}
	stackPars['scale'] = 'normalize_median'
	stackPars['reject'] = 'sigma_clip'
	stackPars['nsplit'] = 8
	stackFun = bokutil.ClippedMeanStack(**stackPars)
	with open(os.path.join('config','illum_list.txt')) as f:
		inputFiles = f.read().splitlines()
	print 'input: ',inputFiles
	print 'output: ',skyFlatFile
	stackFun.stack(inputFiles,skyFlatFile)

# argh. spline requires monotonically increasing coordinates
# as input, but my origin is at the center...
def form_spline_im(ccdName,splinefun,xx,yy):
	if type(ccdName) is int:
		ccdName = 'CCD%d' % (ccdName+1)
	if ccdName == 'CCD1':
		ccdim = splinefun(xx[0,:],yy[:,0]).T
	elif ccdName == 'CCD2':
		ccdim = splinefun(xx[0,:],yy[::-1,0])[:,::-1].T
	elif ccdName == 'CCD3':
		ccdim = splinefun(xx[0,::-1],yy[:,0])[::-1,:].T
	elif ccdName == 'CCD4':
		ccdim = splinefun(xx[0,::-1],yy[::-1,0])[::-1,::-1].T
	return ccdim

def fit_illumination(rmDataMap,nbin=16,asPoly=False,order=3):
	fits = bokutil.BokMefImage(rmDataMap('IllumCorrFlat'),
	                           mask_file=rmDataMap('MasterBadPixMask4'),
	                           read_only=True)
	im = fits.make_fov_image(nbin,'sky',binfunc=np.ma.mean,
	                         binclip=True,single=True)
	if asPoly:
		# fit a polynomial to the binned mosaic image
		poly_model = models.Polynomial2D(degree=order)
		fitfun = fitting.LinearLSQFitter()
		illum = fitfun(poly_model,im['x'],im['y'],im['im'])
	else:
		nIter = 1
		mask = im['im'].mask
		knotints = np.linspace(0.1,1.0,5)
		for iterNum in range(nIter):
			ii = np.where(~mask)
			xmin,xmax = im['x'].min(),im['x'].max()
			tx = np.concatenate([ knotints[::-1]*xmin, knotints*xmax])
			ymin,ymax = im['y'].min(),im['y'].max()
			ty = np.concatenate([ knotints[::-1]*ymin, knotints*ymax])
			illum = LSQBivariateSpline(im['x'][ii],im['y'][ii],
			                           im['im'].data[ii],tx,ty,
			                           kx=order,ky=order)#,eps=1e-8)
			# this forms the spline image for comparison during iteration,
			# not used
			#evalIm = np.dstack([ form_spline_im(ccdNum,illum,
			#                                    im['x'][:,:,ccdNum],
			#                                    im['y'][:,:,ccdNum])
			#                       for ccdNum in range(4) ])
	return illum

def test(skyflatf,bpmaskf,asPoly=False,nbin=16):
	'''Generate an illumination correction image from a sky flat and a
	   bad pixel mask.'''
	import pickle
	class dmap(object):
		def __call__(self,k):
			if k=='IllumCorrFlat':
				return skyflatf
			elif k=='MasterBadPixMask4':
				return bpmaskf
	dm = dmap()
	illum = fit_illumination(dm,nbin=nbin,asPoly=asPoly)
	fits = bokutil.BokMefImage(dm('IllumCorrFlat'),
	                           output_file='testillum.fits',
	                           clobber=True)
	for extName,im,hdr in fits:
		xx,yy = fits.get_xy(extName,'sky')
		if asPoly:
			ccdim = illum(xx,yy)
		else:
			ccdim = form_spline_im(extName,illum,xx,yy)
		ccdim /= float(illum(0,0))
		fits.update(ccdim,hdr)
	fits.close()
	with open("testillum.pkl","wb") as handle:
		pickle.dump(illum,handle)

