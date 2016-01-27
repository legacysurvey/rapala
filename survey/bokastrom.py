#!/usr/bin/env python

import os
import subprocess
import numpy as np
import fitsio
import bass

from astropy import units as u
from astropy.coordinates import SkyCoord

pxscl = 0.45
cfgpath = '/project/projectdirs/desi/users/imcgreer/astrometry-index-4200/cfg'

def solve_bass_image(imagepath,ra=None,dec=None,redo=False,center=False,
                     convert2pv=False,extns='all',writehead=False):
	if not os.path.exists(imagepath):
		print imagepath,' does not exist; skipping'
		return
	if ra is None:
		h = fitsio.read_header(imagepath,0)
		c = SkyCoord(h['RA']+' '+h['DEC'],unit=(u.hourangle,u.deg))
		ra = c.ra.value
		dec = c.dec.value
	sip2pv = os.path.join(os.environ['HOME'],'soft','sip2pv')
	solveargs = ['solve-field','--config',cfgpath]
	solveargs += ['--ra','%.7f'%ra,'--dec','%.7f'%dec,'--radius','1.0']
	solveargs += ['-L','%.2s' % (0.9*pxscl)]
	solveargs += ['-U','%.2s' % (1.1*pxscl)]
	solveargs += ['-u','arcsecperpix']
	solveargs += ['--continue','--no-plots','--no-remove-lines',
	              '--uniformize','0','--no-fits2fits']
	solveargs += ['-N','none','-U','none','-S','none','-M','none',
	              '--rdls','none','--corr','none']
	if center:
		solveargs += ['--crpix-center']
	if convert2pv:
		# in order to do SIP->PV conversion for sextractor;
		# see https://groups.google.com/forum/#!topic/astrometry/qrNW5KtsNCc
		solveargs += ['--tweak-order','4']
	if extns=='all':
		extns = [1,2,3,4]
	for extNum in extns:
		wcsfn = imagepath.replace('.fits','_ext%s.wcs'%extNum)
		solvecmd = solveargs + ['--extension','%d'%extNum]
		if not redo and os.path.exists(wcsfn):
			print wcsfn,' exists, skipping'
			continue
		solvecmd += ['--wcs',wcsfn,imagepath]
		subprocess.call(solvecmd)
		if writehead:
			subprocess.call(['missfits','-c','config/missfits.cfg',imagepath])
		if convert2pv:
			pvimpath = imagepath.replace('.fits','_pv.fits')
			if os.path.exists(pvimpath):
				os.unlink(pvimpath)
			subprocess.call([sip2pv,'-i',imagepath,'-o',pvimpath])

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
#	args = sys.argv[1:]
#	if args[0]=='raw':
#		raw = True
#		args = args[1:]
#	else:
#		raw = False
#	tileId,ditherId = [int(a) for a in args]
#	solve_bass_tile_byid(tileId,ditherId,raw=raw)
	solve_bass_image(sys.argv[1],extns=[1])

