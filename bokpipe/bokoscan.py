#!/usr/bin/env python

import os
import re
from collections import OrderedDict
import numpy as np
from scipy.ndimage.filters import median_filter
from scipy.interpolate import LSQUnivariateSpline
import fitsio

from .bokutil import BokProcess,array_clip,mask_saturation

# argh
ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]
bok90mef_extensions = ['IM%d' % a for a in ampOrder]

def _convertfitsreg(regstr):
	regpattern = r'\[(\d+):(\d+),(\d+):(\d+)\]'
	rv =  [ int(d) for d in  re.match(regpattern,regstr).groups() ]
	# FITS region indices are 1-indexed
	rv[0] -= 1
	rv[2] -= 1
	return rv

def extract_overscan(data,hdr):
	'''Given a 90Prime FITS HDU corresponding to a single amplifier, 
	   extract the overscan regions and trim the image.
	   Returns data, overscan_cols, overscan_rows
	   Output is converted to floats
	'''
	x1,x2,y1,y2 = _convertfitsreg(hdr['BIASSEC'])
	overscan_cols = data[y1:y2,x1:x2].astype(np.float32)
	x1,x2,y1,y2 = _convertfitsreg(hdr['DATASEC'])
	if hdr['NAXIS2'] > y2+1:
		# added overscan rows are not identified in header keywords, just
		# extract any extra rows outside of DATASEC
		overscan_rows = data[y2:,:].astype(np.float32)
	else:
		overscan_rows = None
	data = data[y1:y2,x1:x2].astype(np.float32)
	return ( data,overscan_cols,overscan_rows )

oscan_fit_keywords = ['reject','method','apply_filter','filter_window',
                      'mask_along','clip_iters','clip_sig',
                      'spline_nknots','spline_niter']

def fit_overscan(overscan,**kwargs):
	reject = kwargs.get('reject','sigma_clip')
	method = kwargs.get('method','mean')
	applyFilter = kwargs.get('apply_filter','median')
	windowSize = kwargs.get('filter_window',31)
	maskAlong = kwargs.get('mask_along',[0,1,2,-1])
	along = kwargs.get('along','columns')
	spline_nknots = kwargs.get('spline_nknots',7)
	spline_niter = kwargs.get('spline_niter',2)
	if along == 'rows':
		# make it look like a column overscan for simplicity
		overscan = overscan.transpose()
	npix = overscan.shape[0]
	#
	overscan = np.ma.masked_array(overscan)
	overscan[:,maskAlong] = np.ma.masked
	if along == 'rows':
		# really the columns that need to be masked, but doesn't hurt to
		# mask the edge rows above as well
		overscan[maskAlong,:] = np.ma.masked
	# reject outlier values perpendicular to overscan axis
	if reject == 'sigma_clip':
		overscan = array_clip(overscan,axis=1,**kwargs)
	elif reject == 'minmax':
		overscan[:,overscan.argmax(axis=1)] = np.ma.masked
		overscan[:,overscan.argmin(axis=1)] = np.ma.masked
	# identify outliers along overscan vector (i.e., due to bleed trails
	#   of saturated stars)
	if True:
		overscan = array_clip(overscan,axis=0,
		                      clip_iters=None,clip_sig=3.5,
		                      clip_cenfunc=np.ma.median)
	# collapse overscan into scalar or vector
	if method == 'mean':
		oscan_fit = overscan.mean(axis=1)
	elif method == 'mean_value':
		oscan_fit = np.repeat(overscan.mean(),npix)
	elif method == 'median_value':
		oscan_fit = np.repeat(np.ma.median(overscan),npix)
	elif method == 'cubic_spline':
		# LSQUnivariateSpline only wants the interior knots
		knots = np.linspace(0,npix,spline_nknots)[1:-1]
		oscanvec = overscan.mean(axis=1)
		x = np.arange(npix)
		for iternum in range(spline_niter):
			g = ~oscanvec.mask
			spl_fit = LSQUnivariateSpline(x[g],oscanvec.data[g],t=knots)
			oscan_fit = spl_fit(x)
			res = array_clip(oscanvec-oscan_fit)
			if np.all(res.mask==oscanvec.mask):
				break
			oscanvec.mask |= res.mask
	else:
		raise ValueError
	if 'value' not in method and 'spline' not in method:
		# another round of rejection, this time along the vector
		oscan_fit = array_clip(oscan_fit,axis=0,
		                       clip_iters=None,clip_sig=3.0,
		                       clip_cenfunc=np.ma.mean)
		# fill any vector elements that were rejected as outliers
		oscan_fit = oscan_fit.filled(np.ma.median(oscan_fit))
		# smoothing filter along overscan vector
		if applyFilter == 'median':
			oscan_fit = median_filter(oscan_fit,windowSize)
	return oscan_fit

def overscan_subtract(data,hdr,returnFull=False,row_kwargs=None,**kwargs):
	data,oscan_cols,oscan_rows = extract_overscan(data,hdr)
	colbias = fit_overscan(oscan_cols,**kwargs)
	data[:] -= colbias[:,np.newaxis]
	if oscan_rows is not None:
		# first fit and then subtract the overscan columns at the
		# end of the strip of overscan rows
		# XXX hardcoded geometry
		_colbias = fit_overscan(oscan_rows[:,-20:],**kwargs)
		oscan_rows = oscan_rows[:,:-20] - _colbias[:,np.newaxis]
		# now fit and subtract the overscan rows
		if row_kwargs is None: 
			row_kwargs = {'method':'cubic_spline'}
		rowbias = fit_overscan(oscan_rows,along='rows',**row_kwargs)
		data[:] -= rowbias[np.newaxis,:data.shape[1]]
	else:
		rowbias = None
	if returnFull:
		return data,oscan_cols,oscan_rows,colbias,rowbias
	else:
		return data

class OverscanCollection(object):
	def __init__(self,oscanImgFile,along='columns'):
		self.along = along
		self.imgFile = oscanImgFile
		if along=='columns':
			self.arr_stack = np.hstack
		else:
			self.arr_stack = np.vstack
		self.tmpfn1 = oscanImgFile+'_oscantmp.npy'
		self.tmpfn2 = oscanImgFile+'_restmp.npy'
		if os.path.exists(self.tmpfn1):
			os.unlink(self.tmpfn1)
		if os.path.exists(self.tmpfn2):
			os.unlink(self.tmpfn2)
		self.tmpOscanImgFile = open(self.tmpfn1,'ab')
		self.tmpOscanResImgFile = open(self.tmpfn2,'ab')
		self.files = []
	def close(self):
		os.unlink(self.tmpfn1)
		os.unlink(self.tmpfn2)
	def append(self,oscan,oscanFit,fileName):
		if self.along=='columns':
			resim = (oscan - oscanFit[:,np.newaxis]).astype(np.float32)
		else:
			resim = (oscan - oscanFit[np.newaxis,:]).astype(np.float32)
		try:
			np.save(self.tmpOscanImgFile,oscan.filled(np.nan))
		except:
			np.save(self.tmpOscanImgFile,oscan)
		np.save(self.tmpOscanResImgFile,resim)
		self.files.append(os.path.basename(fileName))
	def write_image(self):
		self.tmpOscanImgFile.close() # could just reset file pointer?
		self.tmpOscanResImgFile.close()
		nfiles = len(self.files)
		if nfiles==0:
			return
		hdr = OrderedDict()
		hdr['NOVSCAN'] = nfiles
		for n,f in enumerate(self.files,start=1):
			hdr['OVSCN%03d'%n] = f
		f1 = open(self.tmpfn1,'rb')
		oscanImg = self.arr_stack([np.load(f1) for i in range(nfiles)])
		f1.close()
		f2 = open(self.tmpfn2,'rb')
		oscanResImg = self.arr_stack([np.load(f2) for i in range(nfiles)])
		f2.close()
		if os.path.exists(self.imgFile+'.fits'):
			os.unlink(self.imgFile+'.fits')
		oscanFits = fitsio.FITS(self.imgFile+'.fits','rw')
		oscanFits.write(oscanImg,header=hdr)
		oscanFits.write(oscanResImg,clobber=False)
		oscanFits.close()
	def n_images(self):
		return len(self.files)

class BokOverscanSubtract(BokProcess):
	_procMsg = 'overscan subtracting %s'
	def __init__(self,**kwargs):
		kwargs.setdefault('header_key','OSCNSUB')
		super(BokOverscanSubtract,self).__init__(**kwargs)
		self.fit_kwargs = { k:v 
		                      for k,v in kwargs.items()
		                        if k in oscan_fit_keywords }
		rkwargs = { k.lstrip('row_'):v 
		              for k,v in kwargs.items()
		                if k.startswith('row_') 
		                  and k.lstrip('row_') in oscan_fit_keywords }
		self.row_fit_kwargs = rkwargs if len(rkwargs)>0 else None
		self.writeOscanImg = kwargs.get('write_overscan_image',False)
		self.oscanColsImgFile = kwargs.get('oscan_cols_file')
		self.oscanRowsImgFile = kwargs.get('oscan_rows_file')
		self.curFileName = None
		self._init_oscan_images()
	def _init_oscan_images(self):
		if self.writeOscanImg:
			self.colImg = { extn:OverscanCollection(self.oscanColsImgFile+
			                                        '_'+extn)
			                   for extn in bok90mef_extensions }
			self.rowImg = { extn:OverscanCollection(self.oscanRowsImgFile+
			                                        '_'+extn,along='rows')
			                   for extn in bok90mef_extensions }
	def _save_oscan_data(self,oscan_cols,colbias,oscan_rows,rowbias,f,extn):
		if self.writeOscanImg:
			self.colImg[extn].append(oscan_cols,colbias,f)
			if oscan_rows is not None:
				self.rowImg[extn].append(oscan_rows,rowbias,f)
	def _finish_oscan_images(self):
		if self.writeOscanImg:
			for extn in bok90mef_extensions:
				self.colImg[extn].write_image()
				self.colImg[extn].close()
				if self.rowImg[extn].n_images() > 0:
					self.rowImg[extn].write_image()
				self.rowImg[extn].close()
	def _preprocess(self,fits,f):
		super(BokOverscanSubtract,self)._preprocess(fits,f)
		self.curFileName = fits.fileName
	def process_hdu(self,extName,data,hdr):
		data,oscan_cols,oscan_rows,colbias,rowbias = \
		         overscan_subtract(data,hdr,returnFull=True,
		                           row_kwargs=self.row_fit_kwargs,
		                           **self.fit_kwargs)
		# write the output file
		hdr['OSCANSUB'] = 'method=%s' % self.fit_kwargs.get('method','default')
		# something changed about what median returns...
		try:
			hdr['OSCANMED'] = float(np.ma.median(colbias).filled(-999))
			if rowbias is not None:
				hdr['OSCNRMED'] = float(np.ma.median(rowbias).filled(-999))
		except:
			hdr['OSCANMED'] = float(np.ma.median(colbias))
			if rowbias is not None:
				hdr['OSCNRMED'] = float(np.ma.median(rowbias))
		self._save_oscan_data(oscan_cols,colbias,oscan_rows,rowbias,
		                      self.curFileName,extName)
		return data,hdr
	def _finish(self):
		self._finish_oscan_images()

class BokOverscanSubtractWithSatFix(BokOverscanSubtract):
	def process_hdu(self,extName,data,hdr):
		data,mask = mask_saturation(extName,data)
		return super(BokOverscanSubtractWithSatFix,self).process_hdu(extName,
		                                                             data,hdr)

