#!/usr/bin/env python

import os,sys
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip

import matplotlib.pyplot as plt

try:
	datadir = os.environ['BASSDIR']
except:
	print 'must set BASSDIR to location of imaging data'
	sys.exit(1)

def fastreadout_analysis():
	getims = lambda f1,f2: [os.path.join(datadir,'20150205',
	                                     'd7058.%04d.fits.gz'%f)
	                          for f in range(f1,f2+1)]
	#
	bias_std_fns = getims(1,10)
	bias_oscan_fns = getims(11,20)
	bias_fast_fns = getims(21,30)
	flat_fastlo_fns = getims(31,40)
	flat_fasthi_fns = getims(41,50)
	flat_stdhi_fns = getims(51,60)
	flat_stdlo_fns = getims(61,62)
	x1,x2 = 300,-300
	y1,y2 = 300,-300
	#
	det = {}
	for readmode in ['standard','fast']:
		det[readmode] = {
		  'gain':np.zeros((16,9)),
		  'readnoise':np.zeros((16,9)),
		}
		if readmode == 'standard':
			biases = bias_std_fns
			flats = flat_stdhi_fns
		else:
			biases = bias_fast_fns
			flats = flat_fasthi_fns
		j = 0
		for b1fn,b2fn,f1fn,f2fn in zip(biases[:-1],biases[1:],
		                               flats[:-1],flats[1:]):
			bias1fits = fits.open(b1fn)
			bias2fits = fits.open(b2fn)
			flat1fits = fits.open(f1fn)
			flat2fits = fits.open(f2fn)
			for ext in range(1,17):
				bias1 = sigma_clip(bias1fits[ext].data[y1:y2,x1:x2])
				bias2 = sigma_clip(bias2fits[ext].data[y1:y2,x1:x2])
				flat1 = sigma_clip(flat1fits[ext].data[y1:y2,x1:x2])
				flat2 = sigma_clip(flat2fits[ext].data[y1:y2,x1:x2])
				_B1 = bias1.mean()
				_B2 = bias2.mean()
				_F1 = flat1.mean()
				_F2 = flat2.mean()
				varF1F2 = (flat1-flat2).var()
				varB1B2 = (bias1-bias2).var()
				# equations from end of sec 4.3 (pg 73) of Howell 2006 
				gain = ( (_F1 + _F2) - (_B1 + _B2) )  / (varF1F2 - varB1B2)
				rdnoise = gain * np.sqrt(varB1B2/2)
				det[readmode]['gain'][ext-1,j] = gain
				det[readmode]['readnoise'][ext-1,j] = rdnoise
				print readmode,j,ext,gain,rdnoise,_F1
			j += 1
	return det

def dump_fastmode_analysis(det):
	print ' '*10,'%6s %6s %6s    ' % ('','gain',''),
	print '%6s %6s %6s  ' % ('','rdnoise','')
	print ' '*10,'%6s %6s %6s   ' % ('std','fast','ratio'),
	print '%6s %6s %6s' % ('std','fast','ratio')
	for i in range(16):
		print 'chip #%2d  ' % (i+1),
		for p in ['gain','readnoise']:
			v1 = sigma_clip(det['standard'][p][i])
			v2 = sigma_clip(det['fast'][p][i])
			print '%6.2f %6.2f %6.3f   ' % \
			       (v1.mean(),v2.mean(),v2.mean()/v1.mean()),
		print

def plot_fastmode_analysis(det):
	for p in ['gain','readnoise']:
		plt.figure(figsize=(7.5,9.5))
		plt.subplots_adjust(0.03,0.03,0.98,0.98)
		for i in range(16):
			plt.subplot(4,4,i+1)
			for mode in ['standard','fast']:
				v = det[mode][p][i]
				eta = {'gain':0.01,'readnoise':0.05}[p]
				bins = np.arange(v.min()-eta,v.max()+2*eta,eta)
				plt.hist(v,bins,histtype='step')

