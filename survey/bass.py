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

tiledb_file = 'bass-newtiles-indesi.fits'
obsdb_file = 'bass-newtiles-observed.fits'

def build_obsdb():
	import glob,re
	dtype = [('utDate','S8'),('fileName','S10'),('expTime','f4'),
	         ('tileId','i4'),('ditherId','i2'),
	         ('ra','f8'),('dec','f8'),('filter','S1')]
	obsdb = np.empty(2e5,dtype=dtype)
	obsfiles = glob.glob(os.path.join(bass_dir,'database','obsed-[gr]-*.txt'))
	nobs = 0
	for fn in sorted(obsfiles):
		_fn = os.path.basename(fn)
		print _fn
		m = re.match(r'obsed-(\w)-(\d+)-(\d+)-(\d+).txt',_fn)
		filt,yr,mo,day = m.groups()
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
					except:
						print 'WARNING: unrecognized tile id %s; skipping' % \
						         tile_s
						continue
					tnum1 = tnum // 10
					tnum2 = tnum % 10
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

def region_tiles(ra1,ra2,dec1,dec2,observed=True):
	if observed:
		tiledb = fits.getdata(os.path.join(bass_dir,obsdb_file))
		ii = np.where((tiledb['ra']>ra1) & (tiledb['ra']<ra2) &
		              (tiledb['dec']>dec1) & (tiledb['dec']<dec2))[0]
	else:
		tiledb = fits.getdata(os.path.join(bass_dir,tiledb_file))
		ii = np.where((tiledb['TRA']>ra1) & (tiledb['TRA']<ra2) &
		              (tiledb['TDEC']>dec1) & (tiledb['TDEC']<dec2))[0]
	return tiledb[ii]

def cfhtw3_tiles(observed=True):
	w3west,w3east = 15*(13.+50/60.), 15*(14+45./60)
	w3south,w3north = 50.7, 56.2
	return region_tiles(w3west,w3east,w3south,w3north,observed=observed)

def ndwfs_tiles(observed=True):
	ndwest,ndeast = 15*14.37, 15*14.62
	ndsouth,ndnorth = 32.5, 36.1
	return region_tiles(ndwest,ndeast,ndsouth,ndnorth,observed=observed)

def check_fields_list():
	files = [ t['utDate']+'/'+t['fileName']+'.fits.gz'
	                 for tiles in [cfhtw3_tiles(),ndwfs_tiles()] 
	                      for t in tiles ]
	with open('checkfields_tiles.txt','w') as f:
		f.write('\n'.join(sorted(files)))


