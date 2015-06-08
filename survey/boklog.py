#!/usr/bin/env python

import os,sys
import glob
import numpy as np
#import fitsio
from astropy.io import fits
from collections import OrderedDict

data_dir = os.environ.get('BOKDATADIR',
                     '/global/project/projectdirs/cosmo/staging/bok')

def boklog(utdate):
	files = glob.glob(data_dir+utdate+'/'+'*.fits.gz')
	files.sort()
	logfn = logdir+'bass.%s.log' % utdate
	if os.path.exists(logfn) or len(files)==0:
		return
	logf = open(logfn,'w')
	nfiles = len(files)
	for i,f in enumerate(files):
		try:
			#h = fitsio.read_header(f,0)
			h = fits.getheader(f)
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
		logf.write('%20s %s %6s %5.1f %7s %15s %6.3f\n' %
		           (fn,h['UT'],filt,h['EXPTIME'],
		            imageType,objname,airmass))
		sys.stdout.write("\r%d/%d" % (i+1,nfiles))
		sys.stdout.flush()
	logf.close()
	print

def log_all():
	utdates = [os.path.basename(d) for d in glob.glob(data_dir+'201?????')]
	utdates.sort()
	for utdate in utdates:
		print 'generating log for ',utdate,'...'
		boklog(utdate)

def load_Bok_log(utdate,logdir='./'):
	'''Load an observing log from a single night.'''
	return np.loadtxt(logdir+'bass.%s.log'%utdate,
	                  dtype=[('fileName','S30'),('utStart','S12'),
	                         ('filter','S8'),('expTime','f4'),('imType','S7'),
	                         ('objectName','S15'),('airmass','f4')])

def load_Bok_logs(logdir='./'):
	'''Load all of the Bok observing logs.'''
	logfiles = sorted(glob.glob(logdir+'*.log'))
	boklogs = OrderedDict()
	for logf in logfiles:
		logfn = os.path.basename(logf)
		utdate = logfn.split('.')[1]
		boklogs[utdate] = load_Bok_log(utdate,logdir)
	return boklogs

if __name__=='__main__':
	log_all()

