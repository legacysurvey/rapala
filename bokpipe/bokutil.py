#!/usr/bin/env python

import os
import fitsio
import numpy as np

def rebin(im,nbin):
	s = np.array(im.shape) / nbin
	return im.reshape(s[0],nbin,s[1],nbin).swapaxes(1,2).reshape(s[0],s[1],-1)

def bok_getxy(hdr,coordsys='image'):
	y,x = np.indices((hdr['NAXIS2'],hdr['NAXIS1']))
	# FITS coordinates are 1-indexed (correct?)
	#x += 1
	#y += 1
	if coordsys == 'image':
		pass
	elif coordsys == 'physical':
		x = hdr['LTM1_1']*(x - hdr['LTV1'])
		y = hdr['LTM2_2']*(y - hdr['LTV2'])
	elif coordsys == 'sky':
		# hacky assumption of orthogonal coordinates but true at this stage
		dx = hdr['CD1_1'] + hdr['CD2_1']
		dy = hdr['CD1_2'] + hdr['CD2_2']
		x = np.sign(dx)*(x - hdr['CRPIX1'])
		y = np.sign(dy)*(y - hdr['CRPIX2'])
	else:
		raise ValueError
	return x,y

class OutputExistsError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class BokMefImage(object):
	'''A wrapper around fitsio that allows the MEF files to be iterated
	   over while updating the data arrays and headers either in-place or
	   to a new file. Also allows for an arbitrary number of masks to be
	   carried with the data.'''
	def __init__(self,fileName,**kwargs):
		self.fileName = fileName
		self.overwrite = kwargs.get('overwrite',False)
		self.keepHeaders = kwargs.get('keep_headers',True)
		self.outFileName = kwargs.get('output_file')
		self.extensions = kwargs.get('extensions')
		self.closeFiles = []
		maskFits = kwargs.get('mask_file')
		if self.outFileName is None:
			self.fits = fitsio.FITS(fileName,'rw')
			self.outFits = self.fits
			self.clobber = True
			self.closeFiles.append(self.fits)
		else:
			self.fits = fitsio.FITS(fileName)
			if os.path.exists(self.outFileName):
				if self.overwrite:
					os.unlink(self.outFileName)
				else:
					raise OutputExistsError
			self.outFits = fitsio.FITS(self.outFileName,'rw')
			hdr = None if not self.keepHeaders else self.fits[0].read_header()
			self.outFits.write(None,header=hdr)
			self.clobber = False
			self.closeFiles.extend([self.fits,self.outFits])
		self.masks = []
		if maskFits is not None:
			self.add_mask(maskFits)
		if self.extensions is None:
			self.extensions = [ h.get_extname() for h in self.fits[1:] ]
		self.curExtName = None
	def add_mask(self,maskFits):
		if type(maskFits) is str:
			maskFits = fitsio.FITS(maskFits)
			self.closeFiles.append(maskFits)
		elif type(maskFits) is not fitsio.fitslib.FITS:
			return ValueError
		self.masks.append(maskFits)
	def update(self,data,header=None):
		self.outFits.write(data,extname=self.curExtName,header=header,
		                   clobber=self.clobber)
	def __iter__(self):
		for self.curExtName in self.extensions:
			data = self.fits[self.curExtName].read()
			hdr = self.fits[self.curExtName].read_header()
			if len(self.masks) > 0:
				mask = self.masks[0][self.curExtName].read().astype(np.bool)
				for m in self.masks[1:]:
					mask |= m[self.curExtName].read().astype(np.bool)
				data = np.ma.masked_array(data,mask=mask)
			yield data,hdr
	def get(self,extName,rows=None,cols=None,header=False):
		if rows is None and cols is None:
			data = self.fits[extName].read()
		else:
			y1,y2 = rows
			x1,x2 = cols
			data = self.fits[extName][y1:y2,x1:x2]
		if header:
			return data,self.fits[extName].read_header()
		else:
			return data
	def close(self):
		for fits in self.closeFiles:
			fits.close()

