#!/usr/bin/env python

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from astropy.stats import sigma_clip
import astropy.io.ascii as ascii_io
import fitsio
import bass

def _convertfitsreg(regstr):
	regpattern = r'\[(\d+):(\d+),(\d+):(\d+)\]'
	rv =  [ int(d) for d in  re.match(regpattern,regstr).groups() ]
	# FITS region indices are 1-indexed
	rv[0] -= 1
	rv[2] -= 1
	return rv

def calc_raw_image_background(imagepath,extNum,margin=800,bmargin=5):
	data,hdr = fitsio.read(imagepath,ext=extNum,header=True)
	x1,x2,y1,y2 = _convertfitsreg(hdr['DATASEC'])
	pix = data[y1+margin:y2-margin,x1+margin:x2-margin].astype(np.float32)
	x1,x2,y1,y2 = _convertfitsreg(hdr['BIASSEC'])
	bias = data[y1+bmargin:y2-bmargin,x1+bmargin:x2-bmargin].astype(np.float32)
	bias = np.median(bias)
	pix -= bias
	pix = sigma_clip(pix,iters=3)
	return dict(medsky=np.ma.median(pix).filled(),meansky=np.ma.mean(pix),
	            rmssky=np.ma.std(pix),bias=bias)

def calc_depth_tile(tile,**kwargs):
	imagepath = os.path.join(bass.bass_data_dir,
	                         tile['utDate'],tile['fileName']+'.fits.gz')
	imstat = calc_raw_image_background(imagepath,
	                                   extNum=kwargs.get('extNum',4),**kwargs)
	return imstat

def calc_depth_all():
	obsdb = bass.load_obsdb()
	outf = open('imagestat.dat','w')
	for ti,tile in enumerate(obsdb):
		print 'tile %d/%d (%s)' % (ti+1,len(obsdb),tile['fileName'])
		if tile['filter']=='g':
			try:
				depthstat = calc_depth_tile(tile)
			except:
				print 'skipping ',tile['fileName']
				continue
			outf.write('%04d %8s %10s %6.1f ' % 
			     (ti,tile['utDate'],tile['fileName'],tile['expTime']))
			outf.write('%8.1f %8.1f %8.1f\n' % 
			     (depthstat['bias'],depthstat['medsky'],depthstat['rmssky']))
	outf.close()

def read_image_stats():
	imstat = ascii_io.read('imagestat.dat',
	                       names=['ti','utDate','fileName','expTime',
	                              'biasADU','skyADU','rmsADU'])
	return imstat

def strip_charts():
	imstat = read_image_stats()
	# XXX use per-image zeropoints once available
	zpt = 25.55
	skymag = -2.5*np.log10(imstat['skyADU']/imstat['expTime']*0.455**-2) + zpt
	plt.figure(figsize=(10,4))
	plt.subplots_adjust(0.06,0.13,0.99,0.98)
	ax1 = plt.subplot(111)
	ax1.plot(skymag)
	utds,ii = np.unique(imstat['utDate'],return_index=True)
	for i,utd in zip(ii,utds):
		offset = 1.5 if np.sum(imstat['utDate']==utd) < 30 else 0.2
		ax1.text(i,skymag[i]-offset,imstat['utDate'][i],size=8,
		         rotation='vertical',ha='left',va='bottom')
	ax1.axvspan(ii[utds==20150305],ii[utds==20150313],color='c',alpha=0.5)
	ax1.axhline(22.1,c='r')
	ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
	ax1.set_ylim(22.5,17)
	ax1.set_xlim(-10,6200)
	ax1.set_xlabel('image number')
	ax1.set_ylabel('g sky brightness [mag/asec2]')

if __name__=='__main__':
	calc_depth_all()

