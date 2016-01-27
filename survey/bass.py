#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
#import astropy.io.ascii as ascii_io

try:
	bass_dir = os.environ['BASSDIR']
	bass_data_dir = os.environ['BASSDATA']
except:
	print 'must set env variables BASSDIR and BASSDATA'
	raise ImportError

try:
	rdxdir = os.environ['BASSRDXDIR']
except:
	try:
		rdxdir = os.path.join(os.environ['GSCRATCH'],'rmreduce')
	except:
		rdxdir = None

tiledb_file = 'bass-newtiles-indesi.fits'
obsdb_file = 'bass-newtiles-observed.fits'

def build_obsdb(noskip=True):
	import glob,re
	dtype = [('utDate','S8'),('fileName','S35'),('expTime','f4'),
	         ('tileId','i4'),('ditherId','i2'),
	         ('ra','f8'),('dec','f8'),('filter','S1')]
	obsdb = np.empty(2e5,dtype=dtype)
	obsfiles_new = glob.glob(os.path.join(bass_dir,'database',
	                                      'obsed-[gr]-*.txt'))
	obsfiles_2015 = glob.glob(os.path.join(bass_dir,'database','2015_old',
	                                       'obsed-[gr]-*.txt'))
	obsfiles = obsfiles_2015 + obsfiles_new
	nobs = 0
	for fn in sorted(obsfiles):
		_fn = os.path.basename(fn)
		print _fn,
		m = re.match(r'obsed-(\w)-(\d+)-(\d+)-(\d+).txt',_fn)
		try:
			filt,yr,mo,day = m.groups()
			print ' loaded'
		except:
			print ' skipped'
			continue
		utd = yr+mo+day
		#obstab = ascii_io.read(fn)
		with open(fn) as f:
			for l in f:
				if l.startswith('#'):
					continue
				bokfn,expt_s,tile_s,ra_s,dec_s = l.strip().split()
				if tile_s.startswith('X'):
					print 'skipping bad tile ',l.strip()
					continue
				elif tile_s.startswith('cal'):
					obsdb['tileId'][nobs] = -99
					obsdb['ditherId'][nobs] = -99
				else:
					try:
						tnum = int(tile_s)
						tnum1 = tnum // 10
						tnum2 = tnum % 10
					except:
						if noskip:
							tnum1,tnum2 = -1,-1
						else:
							print 'WARNING: unrecognized tile id %s; ' \
							      'skipping' % tile_s
					obsdb['tileId'][nobs] = tnum1
					obsdb['ditherId'][nobs] = tnum2
				obsdb['utDate'][nobs] = utd
				m = re.match(r'.*(d\d\d\d\d.\d\d\d\d).*',bokfn).groups()[0]
				obsdb['fileName'][nobs] = m
				obsdb['expTime'][nobs] = float(expt_s)
				obsdb['ra'][nobs] = float(ra_s)
				obsdb['dec'][nobs] = float(dec_s)
				obsdb['filter'][nobs] = filt
				nobs += 1
	obsdb = obsdb[:nobs]
	print 'ingested %d observed tiles' % nobs
	fits.writeto(os.path.join(bass_dir,obsdb_file),obsdb,clobber=True)

def load_tiledb():
	return fits.getdata(os.path.join(bass_dir,tiledb_file))

def load_obsdb():
	return fits.getdata(os.path.join(bass_dir,obsdb_file))

def region_tiles(ra1,ra2,dec1,dec2,observed=True):
	if observed:
		tiledb = load_obsdb()
		ii = np.where((tiledb['ra']>ra1) & (tiledb['ra']<ra2) &
		              (tiledb['dec']>dec1) & (tiledb['dec']<dec2))[0]
	else:
		tiledb = load_tiledb()
		ii = np.where((tiledb['TRA']>ra1) & (tiledb['TRA']<ra2) &
		              (tiledb['TDEC']>dec1) & (tiledb['TDEC']<dec2))[0]
	return tiledb[ii]

def obs_summary(doplot=False,saveplot=False):
	from collections import defaultdict
	tiledb = load_tiledb()
	obsdb = load_obsdb()
	tid = np.array([int(tid) for tid in tiledb['TID']])
	nobs = np.zeros((tiledb.size,3),dtype=int)
	tileList = {1:defaultdict(list),2:defaultdict(list),3:defaultdict(list)}
	for n,row in enumerate(obsdb):
		if row['tileId']>0:
			try:
				i = np.where(tid==row['tileId'])[0][0]
			except:
				print 'tile ',row['tileId'],' is not in db'
				continue
			nobs[i,row['ditherId']-1] += 1
			tileList[row['ditherId']][row['tileId']].append(n)
	print 'total tiles: '
	for i in range(3):
		print 'D%d: %d' % (i+1,np.sum(nobs[:,i]))
	print 'unique tiles: '
	for i in range(3):
		print 'D%d: %d' % (i+1,np.sum(nobs[:,i]>0))
	print 'any dither: ',np.sum(np.any(nobs>0,axis=1))
	print 'repeats: '
	for i in range(3):
		print 'D%d: %d' % (i+1,np.sum(nobs[:,i]>1))
	print 'total repeats: ',np.sum(nobs>1)
	if doplot:
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_pdf import PdfPages
		if saveplot:
			pdf = PdfPages('bass_coverage.pdf')
		for j in range(3):
			fig = plt.figure(figsize=(10,6))
			plt.subplots_adjust(0.03,0.05,0.98,0.95)
			sz = 5 if saveplot else 20
			plt.scatter(tiledb['TRA']/15,tiledb['TDEC'],marker='s',
			            c=np.choose(nobs[:,j],
			                        ['0.9','c','DarkCyan','b','purple','m']),
			            edgecolor='none',s=sz)
			plt.plot([13+5./6,14+45./60,14+45./60,13+5./6,13+5./6],
			         [50.7,50.7,56.2,56.2,50.7],c='k')
			plt.plot([14.37,14.62,14.62,14.37,14.37],
			         [32.5,32.5,36.1,36.1,32.5],c='k')
			plt.xlim(20,5.9)
			plt.ylim(29.7,58)
			plt.title('dither %d total %d unique %d repeats %d' %
			          (j+1,np.sum(nobs[:,j]),np.sum(nobs[:,j]>0),
			           np.sum(nobs[:,j]>1)))
			if saveplot:
				pdf.savefig(fig,orientation='landscape')
		if saveplot:
			pdf.close()
	return nobs,tileList

def nersc_archive_list():
	import fitsio
	from glob import glob
	dirs = sorted(glob(os.path.join(os.environ['BASSDATA'],'BOK_Raw','*')))
	logf = open('nersc_noaoarchive.log','w')
#	errlogf = open('nersc_noaoarchive_errs.log','w')
	for utdir in dirs:
		files = sorted(glob(os.path.join(utdir,'*.fits.fz')))
		print utdir,' %d files' % len(files)
		for f in files:
			h = fitsio.read_header(f,0)
			nersc_path,fn = os.path.split(f)
			nersc_path,nersc_dir = os.path.split(nersc_path)
			try:
				orig_path,orig_fn = os.path.split(h['DTACQNAM'])
			except ValueError:
				orig_path,orig_fn = os.path.split(h['FILENAME'])
#				if not orig_fn.startswith('d7'):
#					errlogf.write('%s missing DTACQNAM/FILENAME\n' % f)
#					continue
			orig_path,orig_dir = os.path.split(orig_path)
			orig_fn = orig_fn.rstrip('.fz')
			exptime = h['EXPTIME']
			imtype = h['IMAGETYP']
			objname = h['OBJECT'].strip()
			if len(objname)==0:
				objname = '<null>'
			logf.write('%8s %30s %18s %10s %6.1f %s\n' %
			           (nersc_dir,fn,orig_fn,imtype,exptime,objname))
		logf.flush()
	logf.close()
#	errlogf.close()

if __name__=='__main__':
	#build_obsdb()
	nersc_archive_list()


