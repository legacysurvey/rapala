#!/usr/bin/env python

import os
import re
import numpy as np
from astropy.stats import sigma_clip
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

if __name__=='__main__':
	calc_depth_all()

