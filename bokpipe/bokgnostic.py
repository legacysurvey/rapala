#!/usr/bin/env python

import os
from collections import defaultdict
from itertools import imap
import numpy as np
import multiprocessing as mp
from functools import partial
import ctypes
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.wcs import WCS
import fitsio

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from . import bokutil
from . import bokproc
from . import bokastrom

def _image_fun(nbin,targets,read_data_fun,f):
	nx,ny,nbinx,nbiny = nbin
	data = read_data_fun(f)
	if data is None:
		return
	xi = (data['x']/nbinx).astype(np.int32)
	yi = (data['y']/nbiny).astype(np.int32)
	ccdi = data['ccdNum'] - 1
	rv = {}
	for t in targets:
		arr = np.zeros((4,ny,nx,3),dtype=np.float32)
		g = np.where(~np.isnan(data[t]))
		v = data[t][g]
		b = np.vstack([v,v**2,np.repeat(1,len(v))]).transpose()
		np.add.at(arr, (ccdi[g],yi[g],xi[g]), b)
		rv[t] = arr
	return rv

class FocalPlaneMap(object):
	def __init__(self,nbin,targets,readDataFun,nproc=1,
	             maxMem=10,prefix='fpmap'):
		try:
			self.nbinx,self.nbiny = nbin
		except:
			self.nbinx = self.nbiny = nbin
		self.nx = 4096 // self.nbinx
		self.ny = 4032 // self.nbiny
		self.N = self.nx * self.ny
		self.targets = targets
		self.nproc = nproc
		if self.nproc == 1:
			self._map = imap
		else:
			self.pool = mp.Pool(self.nproc)
			self._map = self.mp_map
			memsz = self.nx*self.ny*4*3 * len(self.targets) * self.nproc
			memsz = float(memsz*4) / 1024**3
			self.maxChunkSize = int(np.floor(maxMem/memsz))
		self.readDataFun = readDataFun
		self.outIms = None
		self.prefix = prefix
	def mp_map(self,fun,files):
		nChunks = len(files) // (self.nproc*3)
		chunkSize = max(1,min(self.maxChunkSize,nChunks))
		return self.pool.imap_unordered(fun,files,chunkSize)
	def ingest(self,files):
		_fun = partial(_image_fun,(self.nx,self.ny,self.nbinx,self.nbiny),
		               self.targets,self.readDataFun)
		rv = self._map(_fun,files)
		if self.outIms is None:
			self.outIms = rv.next()
		for res in rv:
			if res is None:
				continue
			for t in self.targets:
				self.outIms[t] += res[t]
		for t in self.targets:
			# mean is SUM(x) / N
			self.outIms[t][...,0] /= self.outIms[t][...,2]
			# rms is SQRT( SUM(x^2)/N - <x>^2)
			self.outIms[t][...,1] = np.sqrt( 
			          self.outIms[t][...,1]/self.outIms[t][...,2] - 
			                    self.outIms[t][...,0]**2)
		return self.outIms
	def write(self,prefix=None):
		if not prefix:
			prefix = self.prefix
		for t in self.targets:
			for k,s in enumerate(['mean','rms','n']):
				outfn = prefix+'_%s_%s.fits' % (t,s)
				if os.path.exists(outfn):
					os.remove(outfn)
				fits = fitsio.FITS(outfn,'rw')
				fits.write(None)
				for i in range(4):
					fits.write(self.outIms[t][i,:,:,k],extname='CCD%d'%(i+1))
				fits.close()
	def make_plots(self,prefix=None,**kwargs):
		vmin = kwargs.get('vmin',-0.15)
		vmax = kwargs.get('vmax', 0.15)
		pclip = kwargs.get('pclip',(1,99))
		plt.ioff()
		if not prefix:
			prefix = self.prefix
		pdffn = prefix+'.pdf'
		with PdfPages(pdffn) as pdf:
			for t in self.targets:
				for k,s in enumerate(['mean','rms','n']):
					basenm = prefix+'_%s_%s' % (t,s)
					fitsfn = basenm+'.fits'
					fig = plt.figure(figsize=(7.25,8))
					fig.subplots_adjust(0.02,0.08,0.98,0.94,0.03,0.01)
					cax = fig.add_axes([0.1,0.03,0.8,0.03])
					try:
						fitsim = fitsio.FITS(fitsfn)
					except IOError:
						continue
					ims = np.array([hdu.read() for hdu in fitsim[1:]])
					if False: #s == 'mean':
						_vmin,_vmax = vmin,vmax
					else:
						_vmin,_vmax = np.percentile(ims,pclip)
					for im,ccdnum in zip(ims,[1,3,2,4]): # order as on sky
						ax = fig.add_subplot(2,2,ccdnum)
						_im = ax.imshow(im,vmin=_vmin,vmax=_vmax,
						                cmap=plt.cm.coolwarm,
						                origin='lower',
						                interpolation='nearest')
						ax.xaxis.set_visible(False)
						ax.yaxis.set_visible(False)
						if ccdnum==1:
							cb = fig.colorbar(_im,cax,
							                  orientation='horizontal')
							cb.ax.tick_params(labelsize=9)
					fig.text(0.5,0.98,('%s_%s'%(t,s)).replace('_',' - '),
					         ha='center',va='top',size=14)
					pdf.savefig()
					plt.close()
		plt.ion()

def obs_meta_data(dataMap,outFile='obsmetadata.fits'):
	files_and_frames = dataMap.getFiles(imType='object',with_frames=True)
	cols = defaultdict(list)
	for f,i in zip(*files_and_frames):
		frameId,expTime = dataMap.obsDb['frameIndex','expTime'][i]
		procf = dataMap('sky')(f)
		catf = dataMap('cat')(f)
		try:
			imFits = fitsio.FITS(procf)
			print frameId,procf,' processing'
		except:
			print frameId,procf,' does not exist'
			continue
		cols['frameIndex'].append(frameId)
		hdrs = [ imFits[extName].read_header()
		                 for extName in ['CCD1','CCD2','CCD3','CCD4'] ]
		cols['biasDN'].append([ h['OSCANMED'] for h in hdrs ])
		try:
			cols['skyElPerSec'].append(hdrs[0]['SKYVAL'])
		except:
			cols['skyElPerSec'].append(0)
		cols['avCcdGain'].append([ np.mean([ h['GAIN%02d'%ampNum]
		                          for ampNum in (ccdn*4+np.arange(4)+1) ])
		                            for ccdn,h in enumerate(hdrs) ])
		try:
			cols['fringeScale'].append([ h['FRNGSCL'] for h in hdrs ])
		except:
			cols['fringeScale'].append([0]*4)
		try:
			if 'TPV' in hdrs[0]['CTYPE1']:
				wcshdrs = hdrs
			else:
				aheadf = procf.replace('.fits','.ahead')
				wcshdrs = bokastrom.read_headers(aheadf)
			radec = [ WCS(h).all_pix2world(0,0,1) for h in wcshdrs ]
			radec = np.array(radec)
			raCent,decCent = np.mean(radec,axis=0)
		except:
			print i,procf,' unable to extract WCS'
			raCent,decCent = 0.0,0.0
		cols['raCenter'].append(raCent)
		cols['decCenter'].append(decCent)
		try:
			catFits = fitsio.FITS(catf)
			fwhm = []
			for extNum in range(1,5):
				objs = catFits[extNum].read()
				ii = np.where((objs['FLAGS']==0) &
				              (objs['MAG_AUTO']-objs['MAG_PSF'] > -0.15) &
				              (objs['MAGERR_PSF'] < 0.1))[0]
				if len(ii) > 5:
					v = sigma_clip(objs['FWHM_IMAGE'][ii])
					fwhm.append(float(np.ma.median(v)))
				else:
					fwhm.append(-1)
		except IOError:
			print i,catf,' does not exist'
			fwhm = [-1,-1,-1,-1]
		cols['fwhmPix'].append(fwhm)
	cols = dict(cols)
	tab = Table(cols)
	tab.write(outFile,overwrite=True)

def check_gain_bal(fileName,badPixMaskFile=None,
                   showMean=True,showMode=True,**kwargs):
	maskFits = None
	if badPixMaskFile is not None:
		maskFits = fitsio.FITS(badPixMaskFile)
	gainBalance = bokproc.BokCalcGainBalanceFactors(save_arrays=True,
                                                 mask_map=lambda f: maskFits,
	                                                **kwargs)
	gainBalance.process_files([fileName])
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(16):
		ax = plt.subplot(4,4,i+1)
		extn = 'IM%d'%(i+1)
		pix = gainBalance.arrays[i].flatten()
		modalVal = 3*np.ma.median(pix)-2*pix.mean()
		plt.hist(pix,100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(gainBalance.skyVals[i],color='r',lw=1.2)
		if showMean:
			plt.axvline(pix.mean(),color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,
		         np.ma.median(pix),pix.mean(),
		         pix.std(),pix.min(),pix.max())

def check_sky_level(fileName,maskFile=None,statreg='ccd_central_quadrant',
                    showMean=True,showMode=True,**kwargs):
	fits = fitsio.FITS(fileName)
	maskFits = None
	if maskFile is not None:
		maskFits = fitsio.FITS(maskFile)
	statsPix = bokutil.stats_region(statreg)
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(4):
		ax = plt.subplot(2,2,i+1)
		extn = 'CCD%d'%(i+1)
		pix = fits[extn].read()[statsPix]
		if maskFits is not None:
			pix = np.ma.masked_array(pix,
			                         mask=maskFits[extn].read()[statsPix]>0)
		medVal,rms,pix = bokutil.array_stats(pix,method='median',
		                                     clip=True,rms=True,
		                                     retArray=True,**kwargs)
		meanVal = pix.mean()
		modalVal = 3*medVal-2*meanVal
		plt.hist(pix.flatten(),100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(medVal,color='r',lw=1.2)
		if showMean:
			plt.axvline(meanVal,color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,medVal,meanVal,
		         pix.std(),pix.min(),pix.max())


html_diag_head = '''<html>
<head>
 <style type="text/css">
  body {font-family: monospace}
  body {font-size: xx-small}
  td {padding-right: 10px}
 </style>
<body><table border=1>
'''

html_diag_foot = '''
</table>
</body></html>
'''

def html_table_entry(val,status):
	clrs = {
	  'missing':'black',
	  'nominal':'white',
	  'warning':'yellow',
	  'bad':'red',
	  'weird':'grey',
	}
	clr = clrs.get(status,'white')
	return r'<td bgcolor=%s>%s</td>' % (clr,val)

def run_scamp_diag(imageFiles,ncols=4,**kwargs):
	tabf = open(os.path.join('scamp_diag.html'),'w')
	tabf.write(html_diag_head)
	rowstr = ''
	for n,imFile in enumerate(imageFiles):
		aheadfn = imFile.replace('.fits','.ahead')
		rowstr += r'<td>%s</td>' % os.path.basename(imFile)
		if not os.path.exists(aheadfn):
			print imFile,' missing'
			rowstr += html_table_entry('','missing')
		else:
			hdrs = bokastrom.read_headers(aheadfn)
			rms = [ 3600*np.sqrt(hdr['ASTRRMS1']**2 + hdr['ASTRRMS2']**2)
			            for hdr in hdrs ]
			rms = np.mean(rms)
			if rms < 0:
				status = 'weird'
			elif rms > 0.4:
				status = 'bad'
			elif rms > 0.2:
				status = 'warning'
			else:
				status = 'nominal'
			rowstr += html_table_entry('%.2f'%rms,status)
		if (n%ncols)==ncols-1:
			tabf.write(r'<tr>'+rowstr+r'</tr>'+'\n')
			rowstr = ''
	tabf.write(html_diag_foot)
	tabf.close()

