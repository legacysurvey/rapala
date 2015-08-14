#!/usr/bin/env python

import os,sys
import numpy as np
try:
	import fitsio
except:
	from astropy.io import fits
from astropy.stats import sigma_clip

from ninetyprime import ampOrder,colbias

import matplotlib.pyplot as plt
from matplotlib import ticker

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
	if fn=='bass':
		import basslog
		datadir = get_BASS_datadir()
		logs = basslog.load_Bok_logs('../survey/logs/')
		dpfx = ''
	else:
		import boklog
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

def linearity_check():
	import basslog
	datadir = get_BASS_datadir()
	logs = basslog.load_Bok_logs('../survey/logs/')
	flatlist = [('20150112',1.0),('20150112',2.0),
		        ('20150115',0.6),('20150115',2.0),
		        ('20150205',0.2),('20150205',3.0)]
	logf = open('linearity_check_list.txt','w')
	ims = []
	for utd,texp in flatlist:
		print utd,texp
		ii = np.where((logs[utd]['imType']=='flat') &
		               (np.abs(logs[utd]['expTime']-texp)<0.1) &
		               (logs[utd]['filter']=='g'))[0]
		ii = ii[1:] # skip the first one
		if utd=='20150115' and texp==2.0:
			ii = ii[10:20] # skips over some bad files
		else:
			ii = ii[:10] # limit to 10
		files = []
		for i in ii:
			f = fitsio.FITS(os.path.join(datadir,utd,
		                                  logs[utd]['fileName'][i]+'.fits.gz'))
			files.append(f)
			logf.write('%s %s %.1f\n' % (utd,logs[utd]['fileName'][i],texp))
		data = []
		for ext in range(1,17):
			extdata = [ colbias(f[ext])[0][512:1536,512:1536] for f in files ]
			extdata = np.dstack(extdata)
			data.append(extdata)
		data = np.rollaxis(np.array(data),0,3)
		# this was to look for the alternating shutter speed effect, but it
		# doesn't seem to show up...
		perimage = data.reshape(1024**2,16,-1)
		print np.median(perimage,axis=0)
		# force it to be an even number
		perimage = perimage[...,:(perimage.shape[-1]//2)*2]
		pair_ratio = perimage[...,1::2]/perimage[...,::2]
		pair_ratio = np.median(pair_ratio,axis=0)
		pair_ratio = pair_ratio.mean(axis=0)
		print pair_ratio
		data = np.median(data,axis=-1)
		ims.append(data)
		print
	ims = np.rollaxis(np.array(ims),0,4)
	fratio = ims[...,1::2] / ims[...,::2]
	fratio = fratio.reshape(-1,16,fratio.shape[-1])
	exptimes = np.array([f[1] for f in flatlist])
	expratio = exptimes[1::2]/exptimes[::2]
	outf = open('bok_linearity.dat','w')
	mfratio = np.mean(fratio,axis=0)
	sfratio = np.std(fratio,axis=0)
	for i in range(len(expratio)):
		outf.write('# %.1f %.1f %.1f %.1f\n' % (exptimes[2*i],exptimes[2*i+1],
	                                            ims[...,2*i].mean(),
	                                            ims[...,2*i+1].mean()))
		for j in range(16):
			outf.write('%.3f %.3f\n' % (mfratio[j,i],sfratio[j,i]))
	outf.close()
	logf.close()

def plot_linearity_check():
	h = [l[2:].strip().split() for l in open('bok_linearity.dat').readlines() 
	           if l.startswith('#')]
	h = np.array(h).astype(float)
	expratio = h[:,1]/h[:,0]
	l = np.loadtxt('bok_linearity.dat').reshape(3,16,2)
	plt.figure(figsize=(6,6.5))
	plt.subplots_adjust(0.14,0.09,0.98,0.98,0.25,0.16)
	for i in range(3):
		ax = plt.subplot(3,1,i+1)
		jj = np.argsort(ampOrder)
		plt.errorbar(1+np.arange(16),l[i,jj,0],l[i,jj,1],fmt='bs')
		plt.axhline(np.mean(l[i,:,0]),c='r')
		plt.axhline(expratio[i],c='g')
		ax.text(0.03,0.05,
		        r'$t_1=%.1f,\ \ t_2=%.1f,\ \ counts_1=%d,\ \ counts_2=%d$' %
		        tuple(h[i]),size=11,transform=ax.transAxes)
		ax.text(0.95,0.95,r'counts ratio / exptime ratio = %.3f' %
		        (np.mean(l[i,:,0])/expratio[i]),size=11,ha='right',va='top',
		        transform=ax.transAxes)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
		if i<=1:
			ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
			ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
		else:
			ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
			ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
		plt.xlim(0.1,16.9)
		plt.ylim(0.93*expratio[i],1.03*expratio[i])
	plt.figtext(0.02,0.5,'mean ratio of flats',size=14,
	         ha='left',va='center',rotation='vertical')
	plt.figtext(0.5,0.02,'amplifier number',size=14,ha='center',va='bottom')

def bias_check():
	import basslog
	datadir = get_BASS_datadir()
	logs = basslog.load_Bok_logs('../survey/logs/')
	dtype = [('utDate','S8'),('fileName','S30'),('oscanMean','f4',16),
	         ('meanResidual','f4',16),('rmsResidual','f4',16)]
	biaslog = []
	for utd in sorted(logs.keys()):
		ii = np.where(logs[utd]['imType']=='zero')[0]
		print utd,len(ii)
		for i in ii:
			f = fitsio.FITS(os.path.join(datadir,utd,
		                                  logs[utd]['fileName'][i]+'.fits.gz'))
			biasent = np.zeros(1,dtype=dtype)
			biasent['utDate'] = utd
			biasent['fileName'] = logs[utd]['fileName'][i]
			for ext in range(1,17):
				try:
					im,bias = colbias(f[ext])
				except:
					continue
				bias_residual = sigma_clip(im[512:1536,512:1536])
				biasent['oscanMean'][0,ext-1] = bias
				biasent['meanResidual'][0,ext-1] = bias_residual.mean()
				biasent['rmsResidual'][0,ext-1] = bias_residual.std()
			biaslog.append(biasent)
	biaslog = np.concatenate(biaslog)
	fitsio.write('bass_bias_log.fits',biaslog,clobber=True)

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

def bias_drops():
	import basslog
	from ninetyprime import extract_colbias
	nightlyLogs = basslog.load_Bok_logs('../survey/logs/')
	extNum = 1
	logf = open('bias_drops.txt','w')
	for utd in sorted(nightlyLogs.keys())[1:]:
		print utd
		log = nightlyLogs[utd]
		imdir = os.path.join(os.environ['BASSDATA'],utd)
		for fn,imType in zip(log['fileName'],log['imType']):
			imhdu = fitsio.FITS(os.path.join(imdir,fn+'.fits.gz'))
			im,bias = extract_colbias(imhdu[extNum])
			bias = bias.astype(np.float32)
			centerbias = np.median(bias[500:-500,5:-5])
			bias -= centerbias
			if np.median(bias[5:20,5:-5]) < -15:
				print utd,fn,imType
				logf.write('%s %s %s\n' % (utd,fn,imType))
	logf.close()


if __name__=='__main__':
	#calc_all_gain_rdnoise(10)
	#calc_all_gain_rdnoise(10,'sdss')
	#linearity_check()
	#bias_check()
	bias_drops()

