#!/usr/bin/env python

import os
import re
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack

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

def build_obsdb(noskip=True,update=True,onlygood=True):
	import glob,re,shutil
	from urllib2 import urlopen
	import tarfile
	if update:
		resp = urlopen('http://batc.bao.ac.cn/BASS/lib/exe/fetch.php?media=observation:observation:database.tar.gz')
		f = resp.read()
		tarname = os.path.join(bass_dir,'database.tar.gz')
		outf = open(tarname,'wb')
		outf.write(f)
		outf.close()
		shutil.rmtree(os.path.join(bass_dir,'database'))
		tar = tarfile.open(tarname)
		tar.extractall(path=bass_dir)
		tar.close()
	if onlygood:
		obsfiles = ['obsed-g-2015-good.txt','obsed-g-2016-0102-good.txt',
		            'obsed-r-2015-good.txt','obsed-r-2016-0102-good.txt',]
		obsfiles = [os.path.join(bass_dir,'database',f) for f in obsfiles]
	else:
		obsfiles_new = glob.glob(os.path.join(bass_dir,'database',
		                                      'obsed-[gr]-*.txt'))
		obsfiles_old = glob.glob(os.path.join(bass_dir,'database','201?_old',
		                                      'obsed-[gr]-*.txt'))
		obsfiles = sorted(obsfiles_old + obsfiles_new)
	obsdb = []
	for obsfile in obsfiles:
		# filenames get written in weird ways
		def fnconv(s):
			s1,s2 = re.match('.*\w(\d\d\d\d)[\w.](\d\d\d\d)',s).groups()
			return 'd'+s1+'.'+s2
		# for some reason astropy.Table barfs on reading this in directly
		# so working around it
		def idconv(s):
			try:
				return int(s)
			except:
				return -99
		arr = np.loadtxt(obsfile,dtype=[('fileName','S10'),('expTime','f4'),
		                           ('tileId','i4'),('ra','f8'),('dec','f8')],
		                 converters={0:fnconv,2:idconv})
		if '2015-good' in obsfile:
			# each line in this file is for a single CCD
			arr = arr[::4]
		t = Table(arr)
		t['ditherId'] = t['tileId'] % 10
		t['tileId'] //= 10
		t['filter'] = os.path.basename(obsfile)[6]
		# filename is encoded with last 4 digits of JD
		t['mjd'] = 50000. + np.array([int(d[1:5]) for d in arr['fileName']],
		                             dtype=np.float32)
		obsdb.append(t)
	obsdb = vstack(obsdb)
	obsdb.sort('fileName')
	#
	print 'ingested %d observed tiles' % len(obsdb)
	outf = obsdb_file if onlygood else obsdb_file.replace('.fits','_all.fits')
	obsdb.write(os.path.join(bass_dir,outf),overwrite=True)
	return

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

def obs_summary(filt='g',mjdstart=None,doplot=False,saveplot=False,pltsfx=''):
	from collections import defaultdict
	tiledb = load_tiledb()
	obsdb = load_obsdb()
	print obsdb
	obsdb = obsdb[obsdb['filter']==filt]
	if mjdstart is not None:
		print obsdb['mjd'].min(),obsdb['mjd'].max()
		obsdb = obsdb[obsdb['mjd']>mjdstart]
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
	print 'any pass: ',np.sum(np.any(nobs>0,axis=1))
	print 'repeats: '
	for i in range(3):
		print 'D%d: %d' % (i+1,np.sum(nobs[:,i]>1))
	print 'total repeats: ',np.sum(nobs>1)
	if doplot:
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_pdf import PdfPages
		if saveplot:
			pdf = PdfPages('bass_coverage_%s%s.pdf'%(filt,pltsfx))
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
			plt.title('filter %s pass %d total %d unique %d repeats %d' %
			          (filt,j+1,np.sum(nobs[:,j]),np.sum(nobs[:,j]>0),
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
	build_obsdb()
	#nersc_archive_list()


