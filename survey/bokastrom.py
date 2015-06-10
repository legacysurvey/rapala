#!/usr/bin/env python

import os
import subprocess
import numpy as np
import bass

pxscl = 0.45
cfgpath = '/project/projectdirs/desi/users/imcgreer/astrometry-index-4200/cfg'

rdxdir = os.path.join(os.environ['GSCRATCH'],'rmreduce')

def solve_bass_tile(utDate,fileName,ra,dec,raw=False):
	solveargs = ['solve-field','--config',cfgpath]
	solveargs += ['--ra','%.7f'%ra,'--dec','%.7f'%dec,'--radius','1.0']
	solveargs += ['-L','%.2s' % (0.9*pxscl)]
	solveargs += ['-U','%.2s' % (1.1*pxscl)]
	solveargs += ['-u','arcsecperpix']
	solveargs += ['--continue','--no-plots','--no-remove-lines',
	              '--uniformize','0','--no-fits2fits']
	solveargs += ['--crpix-center','-N','none','-U','none',
	              '-S','none','-M','none',
	              '--rdls','none','--corr','none']
	if raw:
		imagepath = os.path.join(bass.bass_data_dir,
		                         utDate,fileName+'.fits.gz')
		for extNum in [1,5,9,13]:
			wcsfn = 'foo%d.wcs' % extNum
			solvecmd = solveargs + ['--extension','%d'%extNum]
			solvecmd += ['--wcs',wcsfn]
			solvecmd += [imagepath]
			print ' '.join(solvecmd)
	else:
		for ccdNum in range(1,5):
			imagepath = os.path.join(rdxdir,utDate,'ccdproc3',
			                         fileName+'_ccd%d.fits'%ccdNum)
			if not os.path.exists(imagepath):
				print utDate,fileName,' ccd #',ccdNum,
				print ' does not exist; skipping'
				continue
			wcsfn = imagepath.replace('.fits','.wcs')
			if os.path.exists(wcsfn):
				print wcsfn,' exists, skipping'
				continue
			solvecmd = solveargs + ['--wcs',wcsfn,imagepath]
			subprocess.call(solvecmd)
			subprocess.call(['missfits','-c','missfits.cfg',imagepath])

def solve_bass_tile_byid(tileId,ditherId,imageNum=-1,raw=False):
	obsdb = bass.load_obsdb()
	try:
		i = np.where((obsdb['tileId']==tileId) & 
		             (obsdb['ditherId']==ditherId))[0][imageNum]
	except:
		print 'tile %d dither %d has not been observed' % (tileId,ditherId)
		raise ValueError
	return solve_bass_tile(obsdb['utDate'][i],obsdb['fileName'][i],
	                       obsdb['ra'][i],obsdb['dec'][i],raw=raw)

if __name__=='__main__':
	import sys
	args = sys.argv[1:]
	if args[0]=='raw':
		raw = True
		args = args[1:]
	else:
		raw = False
	tileId,ditherId = [int(a) for a in args]
	solve_bass_tile_byid(tileId,ditherId,raw=raw)

