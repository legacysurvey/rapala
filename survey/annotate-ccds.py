#!/usr/bin/env python

from __future__ import print_function

import os
import numpy as np

from astropy.table import Table

def load_etc_results_file(resultsf):
	names = ['fileName','ra','dec','magLimNie','expTime','finalCal',
	         'skyFlux','skyRms','seeing','airmass','filter',
	         'obsCal','obsSkyFlux','obsSkyRms','obsSeeing','obsAirmass',
	         'E(B-V)','calcExpTime',
	         'magLimPoint1','magLimExt1','magLimPoint3','magLimExt3']
	etcVals = Table.read('result.txt',format='ascii',names=names)

def match_results(ccds,etcVals):
	ccds['seeing'][ii] = etcVals['seeing'][jj]
	ccds['zpt'][ii] = etcVals['finalCal'][jj] + 25. \
	                    - 2.5*np.log10(etcVals['expTime'][jj])
	ccds['avsky'][ii] = etcVals['skyFlux'][jj]
	ccds['ebv'][ii] = etcVals['E(B-V)'][jj]

def obslog2ccds(obsLogFile,outfn='bass-ccds-annotated.fits'):
	# read in the observations logs created from trawling the headers,
	# use this as the starting point
	ccds = Table.read(obsLogFile)
	# skip calibration files
	ccds = ccds[ccds['imType']=='object']
	isbass = np.where((ccds['filter']=='g')|(ccds['filter']=='bokr'))[0]
	ccds = ccds[isbass]
	# slice the columns that are needed
	cols = [ 'expTime', 'filter', 'date_obs', 'mjd', 'utObs',
	         'hdrAirmass', 'fileName', 'outsideHumidity', 'outsideTemp', ]
	ccds = ccds[cols]
	# translate some column names to match the DECaLS conventions
	_remap = [ ('expTime','exptime'), ('mjd','mjd_obs'), ('utObs','ut'),
	           ('hdrAirmass','airmass'), ('fileName','image_filename'), 
	           ('outsideHumidity','humidity'), ('outsideTemp','outtemp'),
	  ]
	for m in _remap:
		ccds.rename_column(*m)
	# translate some data values
	ccds['outtemp'] = 5.0/9.0 * (ccds['outtemp']-32) # F->C
	ccds['filter'][ccds['filter']=='bokr'] = 'r'
	# fill some values
	ccds['propid'] = np.repeat('BASS',len(ccds)) # XXX
	ccds['seeing'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['zpt'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['avsky'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['arawgain'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['crpix1'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['crpix2'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['crval1'] = np.zeros(len(ccds),dtype=np.float64)
	ccds['crval2'] = np.zeros(len(ccds),dtype=np.float64)
	ccds['cd1_1'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['cd1_2'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['cd2_1'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['cd2_2'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ebv'] = np.zeros(len(ccds),dtype=np.float32)
	# add dummy values for the per-ccd entries
	ccds['ccdnum'] = np.zeros(len(ccds),dtype=np.int16)
	ccds['ccdname'] = np.repeat('CCDN',len(ccds))
	for amp in ['','a','b','c','d']:
		ccds['ccdzpt'+amp] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ccdphrms'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ccdskyrms'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ccdraoff'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ccddecoff'] = np.zeros(len(ccds),dtype=np.float32)
	ccds['ccdnstar'] = np.zeros(len(ccds),dtype=np.int16)
	for amp in ['','a','b','c','d']:
		ccds['ccdnmatch'+amp] = np.zeros(len(ccds),dtype=np.int16)
	ccds['ccdmdncol'] = np.zeros(len(ccds),dtype=np.float32)
	# now explode the table by 4 copies to make the per-CCD entries
	# write the new table
	ccds.write(outfn,overwrite=True)


if __name__=='__main__':
	import sys
	obslog2ccds(*tuple(sys.argv[1:]))

