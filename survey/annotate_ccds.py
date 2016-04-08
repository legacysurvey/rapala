#!/usr/bin/env python

#from __future__ import print_function

import os
import numpy as np
from astropy.table import Table,join,vstack

from bass import reform_filename,files2tiles,load_obsdb

#def match_results(ccds,etcVals):
#	ccds['seeing'][ii] = etcVals['seeing'][jj]
#	ccds['zpt'][ii] = etcVals['finalCal'][jj] + 25. \
#	                    - 2.5*np.log10(etcVals['expTime'][jj])
#	ccds['avsky'][ii] = etcVals['skyFlux'][jj]
#	ccds['ebv'][ii] = etcVals['E(B-V)'][jj]

def _temp_fn2expid_map(archivelistf):
	fnmap = {}
	with open(archivelistf) as archivef:
		for l in archivef:
			d = l.strip().split()
			noaofn = d[1]#.replace('.fz','').replace('.fits','')
			bokfn = d[2]
			try:
				expid = int(bokfn[1:5]+bokfn[6:10])
			except:
				continue
			fnmap[noaofn] = expid
	return fnmap

def frames2ccds(obsLogFile,outfn='bass-ccds-annotated.fits'):
	# read in the observations logs created from trawling the headers,
	# use the data common to the full frame as the starting point
	frames = Table.read(obsLogFile)
	# identify bass frames
	isbass = ( ((frames['filter']=='g')|(frames['filter']=='bokr')) &
	           (frames['imType']=='object')
	  )
	# XXX need to verify tile id!!!!
	#      this is a hacky way
	isbass &= [ len(s)==5 and s[-1] in '123' for s in frames['objName'] ]
	frames = frames[isbass]
	# slice the columns that are needed
	cols = [ 'expTime', 'filter', 'date_obs', 'mjd', 'utObs',
	         'hdrAirmass', 'fileName', 'outsideHumidity', 'outsideTemp', ]
	frames = frames[cols]
	# translate some column names to match the DECaLS conventions
	_remap = [ ('expTime','exptime'), ('mjd','mjd_obs'), ('utObs','ut'),
	           ('hdrAirmass','airmass'), ('fileName','image_filename'), 
	           ('outsideHumidity','humidity'), ('outsideTemp','outtemp'),
	  ]
	for m in _remap:
		frames.rename_column(*m)
	# translate some data values
	frames['outtemp'] = 5.0/9.0 * (frames['outtemp']-32) # F->C
	frames['filter'][frames['filter']=='bokr'] = 'r'
	frames['propid'] = np.repeat('BASS',len(frames)) # XXX
	frames['image_filename'] = \
	            np.core.defchararray.add(frames['image_filename'],'.fits.fz')
	fnmap = _temp_fn2expid_map('nersc_noaoarchive_thru20160216.log') # XXX
	frames['expnum'] = [ fnmap[fn] for fn in frames['image_filename'] ]
	# allocate dummy entries for per-ccd items
	frames['seeing'] = np.zeros(len(frames),dtype=np.float32)
	frames['zpt'] = np.zeros(len(frames),dtype=np.float32)
	frames['avsky'] = np.zeros(len(frames),dtype=np.float32)
	frames['arawgain'] = np.zeros(len(frames),dtype=np.float32)
	frames['crpix1'] = np.zeros(len(frames),dtype=np.float32)
	frames['crpix2'] = np.zeros(len(frames),dtype=np.float32)
	frames['crval1'] = np.zeros(len(frames),dtype=np.float64)
	frames['crval2'] = np.zeros(len(frames),dtype=np.float64)
	frames['cd1_1'] = np.zeros(len(frames),dtype=np.float32)
	frames['cd1_2'] = np.zeros(len(frames),dtype=np.float32)
	frames['cd2_1'] = np.zeros(len(frames),dtype=np.float32)
	frames['cd2_2'] = np.zeros(len(frames),dtype=np.float32)
	frames['ccdnum'] = np.zeros(len(frames),dtype=np.int16)
	#frames['ccdname'] = np.repeat('CCDN',len(frames))
	frames['ccdphrms'] = np.zeros(len(frames),dtype=np.float32)
	frames['frameskyrms'] = np.zeros(len(frames),dtype=np.float32)
	frames['ccdraoff'] = np.zeros(len(frames),dtype=np.float32)
	frames['ccddecoff'] = np.zeros(len(frames),dtype=np.float32)
	frames['ccdnstar'] = np.zeros(len(frames),dtype=np.int16)
	for amp in ['','a','b','c','d']:
		frames['ccdzpt'+amp] = np.zeros(len(frames),dtype=np.float32)
		frames['ccdnmatch'+amp] = np.zeros(len(frames),dtype=np.int16)
	frames['ccdmdncol'] = np.zeros(len(frames),dtype=np.float32)
	# now explode the framesle by 4 copies to make the per-CCD entries
	ccds = [ frames.copy() for i in range(4) ]
	for j,ccdNum in enumerate(range(1,5)):
		frames['ccdnum'] = ccdNum
	perframekeys = ['seeing',]
	perccdkeys = ['zpt','avsky','crpix1','crpix2','crval1','crval2',
	              'cd1_1','cd1_2','cd2_1','cd2_2','ccdnstar','ccdmdncol',
	              'ccdraoff','ccddecoff']
	perampkeys = ['ccdzpt','ccdnmatch']
	kmap = {'zpt':'zpccd','ccdnstar':'nstarsccd','ccdmdncol':'medgicolor',
	        'ccdraoff':'raoff','ccddecoff':'decoff',
	        'ccdzpt':'zpamp','ccdnmatch':'nstarsamp'}
	for i,frame in enumerate(frames):
		fn = frame['image_filename']
		metaDat = np.load(fn.replace('.fits.fz','.meta.npz'))
		for j,ccdNum in enumerate(range(1,5)):
			# even though these are per-frame still need a copy for each ccd
			for k in perframekeys:
				ccds[j][k][i] = metaDat[k]
			for k in perccdkeys:
				ccds[j][k][i] = metaDat[kmap.get(k,k)][j]
			for ai,aname in enumerate('abcd'):
				for k in perampkeys:
					ccds[j][k+aname][i] = metaDat[kmap.get(k,k)][j,ai]
	allccds = vstack(ccds)
	allccds.sort(['expnum','ccdnum'])
	allccds.write(outfn,overwrite=True)


if __name__=='__main__':
	import sys
	frames2ccds(sys.argv[1])

