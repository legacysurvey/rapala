#!/usr/bin/env python

import os,sys
from glob import glob
import numpy as np
import fitsio
import multiprocessing
from functools import partial
from scipy.ndimage.filters import gaussian_filter
from astropy.stats import sigma_clip
from astropy.table import Table,vstack

from bokpipe.bokoscan import extract_overscan,fit_overscan,overscan_subtract
from bokpipe.bokproc import ampOrder
from bokpipe.bokutil import stats_region,array_clip,array_stats

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
		data = np.empty(1,dtype=[('bias1','S35'),('bias2','S35'),
		                         ('flat1','S35'),('flat2','S35'),
		                         ('biasADU','f4',16),('flatADU','f4',16),
		                      ('biasRmsADU','f4',16),('flatRmsADU','f4',16),
		                         ('gain','f4',16),('rdnoise','f4',16)])
		data['bias1'] = os.path.basename(files[0])
		data['bias2'] = os.path.basename(files[1])
		data['flat1'] = os.path.basename(files[2])
		data['flat2'] = os.path.basename(files[3])
		skip = False
		for ext in range(1,17):
			try:
				bias1,bias2 = [ _data2arr(f[ext].read()[s],100,5000) 
				                  for f in ff[:2]]
				flat1,flat2 = [ _data2arr(f[ext].read()[s],100,50000) 
				                  for f in ff[2:] ]
			except:
				print 'failed with ',ff
				skip = True
				break
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

def find_cal_sequences(log,min_len=5):
	# frameIndex isn't actually a straight table index, but rather a unique
	# id. adding a running index makes the group sorting much clearer.
	t = Table(log).group_by('utDir')
	t['ii'] = np.arange(len(t))
	calseqs = {'zero':[],'flat':[],'zero_and_flat':[]}
	for ut in t.groups:
		iscal = np.where((ut['imType']=='zero')|(ut['imType']=='flat'))[0]
		if len(iscal)==0:
			continue
		# this wouldn't work if someone changed the filter in the middle
		# of a bias sequence... not worth worrying about
		ut_type = ut[iscal].group_by(['imType','filter','expTime'])
		for utt in ut_type.groups:
			if len(utt) < min_len:
				continue
			imType = utt['imType'][0] 
			ii = np.arange(len(utt))
			seqs = np.split(ii,1+np.where(np.diff(utt['frameIndex'])>1)[0])
			seqs = [ np.array(utt['ii'][s]) 
			               for s in seqs if len(s) >= min_len ]
			calseqs[imType].extend(seqs)
	# bias/flat sequences taken in succession, for gain/RN calculation
	# kind of hacky, just look for a set of flats taken roughly close to
	# each set of biases (as in, within 20 minutes)
	max_deltat_minutes = 20.
	bias_times = np.array([ t['mjd'][s[0]] for s in calseqs['zero'] ])
	flat_times = np.array([ t['mjd'][s[0]] for s in calseqs['flat'] ])
	for bt,bs in zip(bias_times,calseqs['zero']):
		j = np.argmin(np.abs(bt-flat_times))
		if 24*60*np.abs(bt-flat_times[j]) < max_deltat_minutes:
			calseqs['zero_and_flat'].append((bs,calseqs['flat'][j]))
	return calseqs

def bias_checks(bias,overscan=False):
	i = 0
	rv = np.zeros(1,dtype=[('fileName','S35'),
	                       ('sliceMeanAdu','f4',(16,)),
	                       ('sliceRmsAdu','f4',(16,)),
	                       ('sliceRangeAdu','f4',(16,)),
	                       ('dropFlag','i4',(16,)),
	                       ('residualMeanAdu','f4',(16,)),
	                       ('residualRmsAdu','f4',(16,))])
	fn = os.path.basename(bias).replace('.fz','').replace('.fits','')
	print 'checking ',fn
	rv['fileName'][i] = fn
	fits = _open_fits(bias)
	for j,hdu in enumerate(fits[1:]):
		imNum = 'IM%d' % ampOrder[j]
		try:
			data = hdu.read().astype(np.float32)
			hdr = hdu.read_header()
		except:
			print 'ERROR: failed to read %s[%d]'%(fn,j+1)
			continue
		if overscan:
			_data,oscan_cols,oscan_rows = extract_overscan(data,hdr)
			#colbias = fit_overscan(oscan_cols,**kwargs)
			if oscan_rows is not None:
				rowbias = fit_overscan(oscan_rows,along='rows',
				                       method='cubic_spline')
				cslice = sigma_clip(data[-22:-2:,2:-22],
				                    iters=1,sigma=3.0,axis=0)
			else:
				cslice = None # No row overscan to check
			bottomslice = data[5:10,-16:-2].mean(axis=0)
			middleslice = data[100:110,-16:-2].mean(axis=0)
		else:
			cslice = sigma_clip(data[1032:1048,2:-22],iters=1,sigma=3.0,axis=0)
			bottomslice = data[5:10,1000:1014].mean(axis=0)
			middleslice = data[100:110,1000:1014].mean(axis=0)
		if cslice is not None:
			cslice = cslice.mean(axis=0)
			cslice = gaussian_filter(cslice,17)
			rv['sliceMeanAdu'][i,j] = cslice.mean()
			rv['sliceRmsAdu'][i,j] = cslice.std()
			rv['sliceRangeAdu'][i,j] = cslice.max() - cslice.min()
		if np.median(middleslice-bottomslice) > 15:
			print 'found drop in ',bias,j
			rv['dropFlag'][i,j] = 1
		bias_residual = overscan_subtract(data,hdr)
		s = stats_region('amp_central_quadrant')
		mn,sd = array_stats(bias_residual[s],method='mean',rms=True,
		                    clip_sig=5.0,clip_iters=2)
		rv['residualMeanAdu'][i,j] = mn
		rv['residualRmsAdu'][i,j] = sd
	return rv

def quick_parallel(fun,input,nproc,**kwargs):
	if nproc > 1:
		fun_with_args = partial(fun,**kwargs)
		p = multiprocessing.Pool(nproc)
		rv = p.map(fun_with_args,input)
		p.close()
		p.join()
	else:
		# instead of creating a 1-process pool just run the sequence
		rv = [ fun(x,**kwargs) for x in input ]
	return np.concatenate(rv)

def run_qa(log,logFits,datadir,nproc=1,dogainrn=True,dobitcheck=True,
	       nsplit=0,nrun=0):
	imType = np.char.rstrip(log['imType'])
	fileNames = np.char.rstrip(log['fileName'])
	utDirs = np.char.rstrip(log['utDir'])
	filePaths = np.char.add(np.char.add(utDirs,'/'),fileNames)
	filePaths = np.char.add(np.char.add(datadir,filePaths),'.fits')
	calseqs = find_cal_sequences(log)
	for imtype in calseqs:
		print 'found %d %s sequences' % (len(calseqs[imtype]),imtype)
	#
	# empirical gain/readnoise calculation
	#
	if dogainrn:
		for bi,fi in calseqs['zero_and_flat']:
			biases = filePaths[bi]
			flats = filePaths[fi]
			# skip the first image in each sequence
			gainrdnoise = calc_gain_rdnoise(biases[1:],flats[1:])
			logFits.write(gainrdnoise,extname='GAINRN')
	#
	# bit integrity check
	#
	if dobitcheck:
		flats = filePaths[np.concatenate(calseqs['flat'])]
		nbits = 8
		bitbit = np.zeros(len(flats),dtype=[('fileName','S35'),
		                                    ('bitFreq','f4',(16,nbits))])
		for i,flat in enumerate(flats):
			fn = os.path.basename(flat).replace('.fz','').replace('.fits','')
			bitbit['fileName'][i] = fn
			fits = _open_fits(flat)
			for j,hdu in enumerate(fits[1:]):
				imNum = 'IM%d' % ampOrder[j]
				data = hdu.read().astype(np.int32)
				npix = float(data.size)
				for bit in range(nbits):
					nbit = np.sum((data&(1<<bit))>0)
					bitbit['bitFreq'][i,j,bit] = nbit/npix
			print 'flat ',i,' out of ',len(flats)
		logFits.write(bitbit,extname='BITCHK')
	#
	# bias ramps
	#
	if len(calseqs['zero'])>0:
		# this checks for bias features using the image region of biases,
		# used to check that the overscan feature search works correctly
		biases = filePaths[np.concatenate(calseqs['zero'])]
		biasrmp = quick_parallel(bias_checks,biases,nproc,overscan=False)
		logFits.write(biasrmp,extname='BIASCHK')
	ii = np.where((imType=='zero')|(imType=='object'))[0]
	if nsplit > 0:
		print 'splits: ',ii[0],len(ii),
		ii = np.array_split(ii,nsplit)[nrun]
		print ii[0],len(ii)
	print 'checking overscans for ',len(ii),' images'
	images = filePaths[ii]
	biasrmp = quick_parallel(bias_checks,images,nproc,overscan=True)
	biasrmp = np.lib.recfunctions.append_fields(biasrmp,'imType',imType[ii],
	                                            dtypes=imType.dtype)
	logFits.write(biasrmp,extname='OSCANCHK')

def run_nightly_checks(utdir,logdir,datadir,redo=False):
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
	run_qa(log,logFits,datadir)
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

def calc_overheads(logdata):
	dt = np.diff(logdata['mjd'])*24*3600 - logdata['expTime'][:-1]
	imt = {'zero':0,'dark':1,'flat':2,'object':3}
	imts = [imt[img['imType'].strip()] for img in logdata[:-1]]
	return imts,logdata['mjd'][:-1],dt

def overhead_report(oheads):
	if type(oheads) is str:
		oheads = np.loadtxt(oheads,unpack=True)
	imt,mjd,dt = oheads
	imt = imt.astype(np.int32)
	plt.scatter(mjd,dt,c=np.choose(imt,['gray','black','cyan','blue']))

def combined_report(logdir,utdates):
	if utdates is None:
		utdates = '*'
	logfiles = sorted(glob(os.path.join(logdir,'log_ut%s.fits'%utdates)))
	print logfiles
	data1,data2 = [],[]
	utbreaks1,utbreaks2 = [0,],[0,]
	oheads = []
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
		log = fits[1].read()
		oheads.append(calc_overheads(log))
	data1 = vstack(data1)
	data2 = vstack(data2)
	oheads = np.hstack(oheads)
	np.savetxt('bass_overheads.txt',oheads.transpose())
	with open('bass_summary.txt','w') as outf:
		gainrn_report(data1,outf,utbreaks1[1:])
		bit_report(data2,outf,utbreaks2[1:])


if __name__=='__main__':
	import argparse
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
	parser.add_argument("-o","--output",type=str,
	                    default="bassqa.fits",
	                    help="output file")
	parser.add_argument("--logfile",type=str,
	                    help="input log file")
	parser.add_argument("--nogainrn",action='store_true',
	                    help="skip gain/readnoise calc (SLOW)")
	parser.add_argument("--nobitcheck",action='store_true',
	                    help="skip bit integrity check")
	parser.add_argument("--nproc",type=int,default=1,
	                    help="set number of processes to run [default 1]")
	parser.add_argument("-R","--redo",action='store_true',
	                    help="ignore existing data and redo")
	parser.add_argument("--numsplit",type=int,default=0,
	                    help="number of chunks to split data into [none]")
	parser.add_argument("--splitnum",type=int,default=0,
	                    help="which chunk number to run")
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
			if args.logfile:
				log = fitsio.read(args.logfile,1)
				if os.path.exists(args.output):
					os.remove(args.output)
				logFits = fitsio.FITS(args.output,'rw')
				logFits.write(None)
				run_qa(log,logFits,args.datadir,nproc=args.nproc,
				       dogainrn=(not args.nogainrn),
				       dobitcheck=(not args.nobitcheck),
				       nsplit=args.numsplit,nrun=args.splitnum)
				logFits.close()
				# if this isn't here multiprocess gets stuck in an infinite
				# loop... why?
				sys.exit(0)
			else:
				run_nightly_checks(utdir,args.logdir,args.datadir,
				                   redo=args.redo)
	if args.report:
		combined_report(args.logdir,args.utdate)

