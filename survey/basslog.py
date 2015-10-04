#!/usr/bin/env python

import os,sys
import glob
import numpy as np
#import fitsio
from astropy.io import fits
from collections import OrderedDict

data_dir = os.environ.get('BOKDATADIR',
                     '/global/project/projectdirs/cosmo/staging/bok')

def boklog(logfn,files):
	logf = open(logfn,'w')
	nfiles = len(files)
	for i,f in enumerate(files):
		try:
			#h = fitsio.read_header(f,0)
			h = fits.getheader(f)
			h1 = fits.getheader(f,1)
		except:
			print 'ERROR: failed to open file: ',f
			continue
		fn = os.path.basename(f).rstrip('.fits.gz')
		imageType = h['IMAGETYP'].strip()
		filt = h['FILTER'].strip()
		objname = h['OBJECT']
		if len(objname.strip())==0:
			objname = 'null'
		# airmass values in header are not very accurate (one decimal place)
		#if not 'AIRMASS' in h:
		#	airmass=-9999
		#else:
		#	airmass=h['AIRMASS']
		if not 'ELEVAT' in h:
			airmass=-9999
		else:
			airmass=1/np.cos(np.radians(90-h['ELEVAT']))
			if np.abs(airmass-h['AIRMASS']) > 0.1:
				print airmass,h['AIRMASS'],h['ELEVAT']
				raise ValueError
		ccdsum = h['CCDSUM'].strip()
		novscanrows = h1['OVRSCAN2']
		logf.write('%10s %14s %6s %5.1f %7s %25s %6.3f %4s %3d\n' %
		           (fn,h['UT'],filt,h['EXPTIME'],
		            imageType,objname,airmass,ccdsum,novscanrows))
		sys.stdout.write("\r%d/%d" % (i+1,nfiles))
		sys.stdout.flush()
	logf.close()
	print

def basslog(utdate,logdir='./'):
	files = glob.glob(os.path.join(data_dir,utdate,'d????.????.fits.gz'))
	files.sort()
	logfn = logdir+'bass.%s.log' % utdate
	if os.path.exists(logfn) or len(files)==0:
		return
	boklog(logfn,files)

def log_all(logdir='./'):
	utdates = [os.path.basename(d) 
	                for d in glob.glob(os.path.join(data_dir,'201?????'))]
	utdates.sort()
	for utdate in utdates:
		print 'generating log for ',utdate,'...'
		basslog(utdate,logdir)

def load_Bok_log(utdate,logdir='./'):
	'''Load an observing log from a single night.'''
	return np.loadtxt(logdir+'bass.%s.log'%utdate,
	                  dtype=[('fileName','S10'),('utStart','S12'),
	                         ('filter','S8'),('expTime','f4'),('imType','S7'),
	                         ('objectName','S25'),('airmass','f4')])

def load_Bok_logs(logdir='./'):
	'''Load all of the Bok observing logs.'''
	logfiles = sorted(glob.glob(logdir+'*.log'))
	boklogs = OrderedDict()
	for logf in logfiles:
		logfn = os.path.basename(logf)
		utdate = logfn.split('.')[1]
		boklogs[utdate] = load_Bok_log(utdate,logdir)
	return boklogs

def overhead_times(utdate,logdir='./'):
	import bass
	from astropy.time import Time
	log = load_Bok_log(utdate,logdir)
	lastt = None
	lastexp = 0
	outf = open('timelog_%s.txt'%utdate,'w')
	for frame in log[:100]:
		h0 = fits.getheader(os.path.join(bass.bass_data_dir,
		                             utdate,frame['fileName']+'.fits.gz'),0)
		ut = h0['UT']
		ut = '-'.join([utdate[:4],utdate[4:6],utdate[6:]])+' '+ut
		t = Time(ut,format='iso',scale='utc')
		if lastt is not None:
			dt = (t-lastt).sec
			outf.write('%s %10s %10.1f %8.1f %8.1f\n' % 
			        (frame['fileName'],frame['imType'],dt,lastexp,dt-lastexp))
		lastt = t
		lastexp = frame['expTime']
	outf.close()

if __name__=='__main__':
	#basslog('20150909','logs/')
	#log_all()
	#for utd in ['20150107','20150117','20150204','20150205','20150211'
	#            '20150305','20150324','20150419','20150426']:
	#	overhead_times(utd,'logs/')
	if True:
		files = glob.glob(os.path.join(data_dir,'20150909','bss','20150908',
		                  'd????.????.fits.gz'))
		files.sort()
		logfn = 'logs/'+'bss.%s.log' % '20150909'
		boklog(logfn,files)

