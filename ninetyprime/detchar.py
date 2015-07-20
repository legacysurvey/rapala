#!/usr/bin/env python

import os,sys
import numpy as np
try:
	import fitsio
except:
	from astropy.io import fits
from astropy.stats import sigma_clip

import matplotlib.pyplot as plt

def _data2arr(data,minVal=0,maxVal=65335):
	return np.ma.masked_array(data.astype(np.float32),
	                          mask=((data<minVal)|(data>maxVal)))

def calc_gain_rdnoise(biases,flats,margin=500):
	rv = []
	x1,x2 = margin,-margin
	y1,y2 = margin,-margin
	for files in zip(biases[:-1],biases[1:],flats[:-1],flats[1:]):
		try:
			ff = [fitsio.FITS(f) for f in files]
		except NameError:
			ff = [fits.open(f) for f in files]
		data = np.empty(1,dtype=[('gain','f4',16),('rdnoise','f4',16)])
		for ext in range(1,17):
			bias1,bias2 = [ _data2arr(f[ext][y1:y2,x1:x2],500,5000) 
			                  for f in ff[:2]]
			flat1,flat2 = [ _data2arr(f[ext][y1:y2,x1:x2],500,40000) 
			                  for f in ff[2:] ]
			_B1 = bias1.mean()
			_B2 = bias2.mean()
			_F1 = flat1.mean()
			_F2 = flat2.mean()
			try:
				varF1F2 = (flat1-flat2).var()
				varB1B2 = (bias1-bias2).var()
			except:
				# some images have wrong format...
				data['gain'][0,ext-1] = -1
				data['rdnoise'][0,ext-1] = -1
				continue
###			print '%s %s %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f' % (os.path.basename(files[0]),os.path.basename(files[1]),ext,_B1,_B2,_F1,_F2,varF1F2,varB1B2,varF1F2-varB1B2)
			# equations from end of sec 4.3 (pg 73) of Howell 2006 
			gain = ( (_F1 + _F2) - (_B1 + _B2) )  / (varF1F2 - varB1B2)
			rdnoise = gain * np.sqrt(varB1B2/2)
			data['gain'][0,ext-1] = gain
			data['rdnoise'][0,ext-1] = rdnoise
		rv.append(data)
	return np.concatenate(rv)


#
# Analysis of data obtained for BASS imaging survey
#

def get_BASS_datadir():
	try:
		datadir = os.environ['BASSDATA']
	except:
		print 'must set BASSDATA to location of imaging data'
		raise ValueError
	return datadir

def calc_all_gain_rdnoise(nmax=5,fn='bass'):
	import boklog
	if fn=='bass'
		datadir = get_BASS_datadir()
		logs = boklog.load_Bok_logs('../survey/logs/')
		dpfx = ''
	else:
		datadir = os.environ['BOK90PRIMERAWDIR']
		logs = boklog.load_Bok_logs()
		dpfx = 'ut'
	utds = sorted(logs.keys())
	detData = []
	fileData = []
	for utd in utds:
		if utd=='20131222': continue # something amiss with these
		# always skip the first three
		ii1 = np.where(logs[utd]['imType']=='zero')[0][3:]
		if len(ii1)>nmax:
			ii1 = ii1[:nmax]
		elif len(ii1)==0:
			continue
		biases = [os.path.join(datadir,dpfx+utd,
		                       logs[utd]['fileName'][i]+'.fits.gz')
		           for i in ii1]
		# always skip the first three
		ii2 = np.where((logs[utd]['imType']=='flat') &
		               (logs[utd]['expTime']>1.0) &
		               (logs[utd]['filter']=='g'))[0][3:]
		if len(ii2)>nmax:
			ii2 = ii2[:nmax]
		elif len(ii2)==0:
			continue
		flats = [os.path.join(datadir,dpfx+utd,
		                      logs[utd]['fileName'][i]+'.fits.gz')
		           for i in ii2]
		if utd == '20150205':
			# hack because this night had mixed readout modes
			getims = lambda f1,f2: [os.path.join(datadir,'20150205',
			                                     'd7058.%04d.fits.gz'%f)
			                          for f in range(f1,f2+1)]
			biases = getims(1,10)
			flats = getims(51,60)
		print 'utd: ',utd,len(ii1),len(ii2)
		c = calc_gain_rdnoise(biases,flats)
		n = c['gain'].shape[0]
		detData.append(c)
		fileData.extend([(utd,biasfn,flatfn) for biasfn,flatfn in
		                        zip(logs[utd]['fileName'][ii1[:n]],
		                            logs[utd]['fileName'][ii2[:n]])])
	detData = np.concatenate(detData)
	fileData = np.array(fileData,dtype=[('utDate','S8'),
	                                 ('biasFn','S10'),('flatFn','S10')])
	fitsio.write('bok90_%s_char.fits'%fn,detData,clobber=True)
	fitsio.write('bok90_%s_char.fits'%fn,fileData)

def fastreadout_analysis():
	datadir = get_BASS_datadir()
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


if __name__=='__main__':
	#calc_all_gain_rdnoise(10)
	calc_all_gain_rdnoise(10,'sdss')

