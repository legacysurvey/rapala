#!/usr/bin/env python

import os,sys
from glob import glob
import numpy as np
import fitsio
from astropy.stats import sigma_clip
from astropy.table import Table,vstack

from bokpipe.bokproc import ampOrder
from bokpipe.bokutil import stats_region,array_clip

import matplotlib.pyplot as plt
from matplotlib import ticker

def _data2arr(data,minVal=0,maxVal=65335,clip=True):
	arr = np.ma.masked_array(data.astype(np.float32),
	                         mask=((data<minVal)|(data>maxVal)))
	return array_clip(arr,clip_iters=2,clip_sig=5.0)

def _open_fits(f):
	try:
		return fitsio.FITS(f)
	except IOError:
		return fitsio.FITS(f+'.fz')

def calc_gain_rdnoise(biases,flats):
	rv = []
	s = stats_region('amp_corner_ccdcenter_1024')
	for files in zip(biases[:-1],biases[1:],flats[:-1],flats[1:]):
		ff = [_open_fits(f) for f in files]
		data = np.empty(1,dtype=[('bias1','S15'),('bias2','S15'),
		                         ('flat1','S15'),('flat2','S15'),
		                         ('biasADU','f4',16),('flatADU','f4',16),
		                      ('biasRmsADU','f4',16),('flatRmsADU','f4',16),
		                         ('gain','f4',16),('rdnoise','f4',16)])
		data['bias1'] = os.path.basename(files[0])
		data['bias2'] = os.path.basename(files[1])
		data['flat1'] = os.path.basename(files[2])
		data['flat2'] = os.path.basename(files[3])
		skip = False
		for ext in range(1,17):
			bias1,bias2 = [ _data2arr(f[ext].read()[s],100,5000) 
			                  for f in ff[:2]]
			flat1,flat2 = [ _data2arr(f[ext].read()[s],100,50000) 
			                  for f in ff[2:] ]
			_B1 = bias1.mean()
			_B2 = bias2.mean()
			_F1 = flat1.mean()
			_F2 = flat2.mean()
			if ext==1 and (np.abs(_F1-_F2)/_F1) > 0.05:
				# a very large jump in the flat field value will throw
				# this calculation off, restrict it to 5% variation
				skip = True
				break
			try:
				varF1F2 = (flat1-flat2).var()
				varB1B2 = (bias1-bias2).var()
			except:
				# some images have wrong format...
				data['gain'][0,ext-1] = -1
				data['rdnoise'][0,ext-1] = -1
				continue
			# equations from end of sec 4.3 (pg 73) of Howell 2006 
			gain = ( (_F1 + _F2) - (_B1 + _B2) )  / (varF1F2 - varB1B2)
			rdnoise = gain * np.sqrt(varB1B2/2)
			data['biasADU'][0,ext-1] = _B1
			data['flatADU'][0,ext-1] = _F1
			data['biasRmsADU'][0,ext-1] = bias1.std()
			data['flatRmsADU'][0,ext-1] = flat1.std()
			data['gain'][0,ext-1] = gain
			data['rdnoise'][0,ext-1] = rdnoise
		if not skip:
			rv.append(data)
	return np.concatenate(rv)


def calc_all_gain_rdnoise(nmax=5,fn='bass'):
	# XXX if keeping then update to use dataMap
	raise NotImplementedError
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

def bias_check():
	# XXX if keeping update to use bokproc
	raise NotImplementedError
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
	# XXX if keeping update to use bokproc
	raise NotImplementedError
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

def nightly_checks(utdir,logdir,redo=False):
	from bokpipe.bokobsdb import generate_log
	print 'running nightly check on ',utdir
	logf = os.path.join(logdir,'log_ut%s.fits' % os.path.basename(utdir))
	if os.path.exists(logf):
		# assume that if log file exists processing is already done
		if not redo:
			return
		else:
			# start the file over
			log = fitsio.read(logf,1)
			os.remove(logf)
			logFits = fitsio.FITS(logf,'rw')
			logFits.write(log)
	else:
		generate_log([utdir],logf)
		logFits = fitsio.FITS(logf,'rw')
		log = logFits[1].read()
	imType = np.char.rstrip(log['imType'])
	fileNames = np.char.rstrip(log['fileName'])
	# find first bias sequence
	biases = []
	is_bias = (imType == 'zero')
	ii = np.where(is_bias)[0]
	for _i,i in enumerate(ii):
		for j in range(_i+1,len(ii)):
			if ii[j] != ii[j-1]+1:
				break
		if (j-_i > 5):
			biases = ii[:j-_i]
			break
	# find first flat sequence in any filter
	flats = []
	is_flat = (imType == 'flat')
	ii = np.where(is_flat)[0]
	for _i,i in enumerate(ii):
		for j in range(_i+1,len(ii)):
			if ( (ii[j] != ii[j-1]+1) or
			     (log['filter'][ii[j]] != log['filter'][ii[j-1]]) or
			     (log['expTime'][ii[j]] != log['expTime'][ii[j-1]]) ):
				break
		if (j-_i > 5):
			flats = ii[:j-_i]
			break
	#
	# empirical gain/readnoise calculation
	#
	biases = [ os.path.join(utdir,f+'.fits') for f in fileNames[biases] ]
	flats = [ os.path.join(utdir,f+'.fits') for f in fileNames[flats] ]
	if len(biases)>3 and len(flats)>3:
		# skip the first image in each sequence
		gainrdnoise = calc_gain_rdnoise(biases[1:],flats[1:])
		logFits.write(gainrdnoise)
	#
	# bit integrity check
	#
	nbits = 8
	bitbit = np.zeros(len(flats),dtype=[('fileName','S15'),
	                                    ('bitFreq','f4',(16,nbits))])
	for i,flat in enumerate(flats):
		bitbit['fileName'][i] = os.path.basename(flat)
		fits = _open_fits(flat)
		for j,hdu in enumerate(fits[1:]):
			imNum = 'IM%d' % ampOrder[j]
			data = hdu.read().astype(np.int32)
			npix = float(data.size)
			for bit in range(nbits):
				nbit = np.sum((data&(1<<bit))>0)
				bitbit['bitFreq'][i,j,bit] = nbit/npix
	logFits.write(bitbit)
	# finish up
	logFits.close()

def _reportfig_init8():
	fig = plt.figure(figsize=(7.5,9.5))
	fig.subplots_adjust(0.1,0.025,0.95,0.975,0.35,0.3)
	return fig

def gainrn_report(data,outf,utbreaks=None,save=True):
	plot_fields = [ ('biasADU',(700,1200), 'biasRmsADU',(2,10)), 
	                ('flatADU',(1e4,6e4), 'flatRmsADU',(0,600)), 
	                ('gain',(0,3), 'rdnoise',(0,15)), 
	]
	for i in range(len(data)):
		outf.write('%4d %s %s %s %s\n' % (i,data['bias1'][i],data['bias2'][i],
		                                  data['flat1'][i],data['flat2'][i]))
	outf.write('\n')
	for fields in plot_fields:
		fig = _reportfig_init8()
		f1,fr1,f2,fr2 = fields
		for i in range(4):
			ax1 = fig.add_subplot(4,2,i+1)
			ax2 = fig.add_subplot(4,2,i+1+4)
			for utb in utbreaks:
				for ax in [ax1,ax2]:
					ax.axvline(utb,c='gray',ls=':',lw=0.5)
			for j in range(4):
				ax1.plot(data[f1][:,4*i+j],label='IM%d'%ampOrder[4*i+j])
				ax2.plot(data[f2][:,4*i+j],label='IM%d'%ampOrder[4*i+j])
			ax1.legend(ncol=4,frameon=False,fontsize=9,columnspacing=1.0)
			ax2.legend(ncol=4,frameon=False,fontsize=9,columnspacing=1.0)
			ax1.set_ylim(*fr1)
			ax2.set_ylim(*fr2)
			ax1.set_ylabel(f1+' CCD%d'%(i+1),size=11)
			ax2.set_ylabel(f2+' CCD%d'%(i+1),size=11)
		for f in [f1,f2]:
			for j in range(4):
				outf.write('%-12s CCD%d\n' % (f,j+1))
				for i in range(len(data)):
					outf.write('%4d  ' % (i))
					for k in range(4):
						outf.write('%8.2f '%data[f][i,4*j+k])
					outf.write('\n')
				outf.write('\n')
		if save:
			plt.savefig('bass_summary_%s.png'%f1[:4])

def _reportfig_init16():
	fig = plt.figure(figsize=(7.5,9.5))
	fig.subplots_adjust(0.1,0.025,0.95,0.975,0.35,0.3)
	return fig

def bit_report(data,outf,utbreaks=None,save=True):
	fig = _reportfig_init16()
	for i in range(16):
		ax = fig.add_subplot(4,4,i+1)
		if utbreaks is not None:
			for utb in utbreaks:
				ax.axvline(utb,c='gray',ls=':',lw=0.5)
		haslab = False # suppresses a warning message
		for j in range(6):
			l = 'bit%d'%j if j//2==i else None
			if l is not None: haslab = True 
			ax.plot(data['bitFreq'][:,i,j],label=l)
		for n in range(data['bitFreq'].shape[0]):
			outf.write('%s  ' % data['fileName'][n])
			outf.write(('%.3f  '*6) % tuple(data['bitFreq'][n,i,:6]))
			outf.write('\n')
		if haslab:
			ax.legend(ncol=2,frameon=False,fontsize=8,columnspacing=1.0)
		ax.set_ylim(0.35,0.65)
		ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
	if save:
		plt.savefig('bass_summary_%s.png'%'bits')

def nightly_report(logf):
	fits = fitsio.FITS(logf)
	data = fits[2].read()
	gainrn_report(data)

def combined_report(logdir):
	logfiles = sorted(glob(os.path.join(logdir,'log_*.fits')))
	data1,data2 = [],[]
	utbreaks1,utbreaks2 = [0,],[0,]
	for i,logf in enumerate(logfiles):
		fits = fitsio.FITS(logf)
		try:
			data1.append(Table(fits[2].read()))
			utbreaks1.append(len(data1[-1])+utbreaks1[-1])
		except:
			pass
		try:
			data2.append(Table(fits[3].read()))
			utbreaks2.append(len(data2[-1])+utbreaks2[-1])
		except:
			pass
	data1 = vstack(data1)
	data2 = vstack(data2)
	with open('bass_summary.txt','w') as outf:
		gainrn_report(data1,outf,utbreaks1[1:])
		bit_report(data2,outf,utbreaks2[1:])


if __name__=='__main__':
	import argparse
	#calc_all_gain_rdnoise(10)
	#calc_all_gain_rdnoise(10,'sdss')
	#linearity_check()
	#bias_check()
	#bias_drops()
	plt.ioff()
	parser = argparse.ArgumentParser()
	parser.add_argument("-n","--nightly",action='store_true',
	                    help="run nightly processing")
	parser.add_argument("-r","--report",action='store_true',
	                    help="make nightly report")
	parser.add_argument("-d","--datadir",type=str,
	                    default="/data/primefocus/bass/",
	                    help="top-level data directory")
	parser.add_argument("-l","--logdir",type=str,
	                    default="/home/mcgreer/basslogs/",
	                    help="log file directory")
	parser.add_argument("-R","--redo",action='store_true',
	                    help="ignore existing data and redo")
	parser.add_argument("-u","--utdate",type=str,
	                    help="restrict UT date")
	args = parser.parse_args()
	#
	if args.utdate is not None:
		utdirs = sorted(glob(os.path.join(args.datadir,args.utdate)))
	else:
		utdirs = sorted(glob(os.path.join(args.datadir,'201?????')))
	if args.nightly:
		for utdir in utdirs:
			nightly_checks(utdir,args.logdir,redo=args.redo)
	if args.report:
		combined_report(args.logdir)

