#!/usr/bin/env python

import os
import shutil
import glob
import subprocess
from copy import copy
import numpy as np
from astropy.io import fits

def scamp_solve(imageFile,catFile,refStarCatFile=None,
                filt='g',savewcs=False,clobber=False,
                check_plots=False,verbose=0,**kwargs):
	configDir = os.path.join(os.path.split(__file__)[0],'config')
	headf = catFile.replace('.fits','.head') 
	wcsFile = imageFile.replace('.fits','.ahead')
	if not clobber and os.path.exists(wcsFile):
		if verbose > 0:
			print wcsFile,' already exists, skipping'
		return
	if refStarCatFile is not None:
		refCatPath = os.path.dirname(refStarCatFile)
		if len(refCatPath)==0:
			refCatPath = '.'
	#
	scamp_cmd_base = ['scamp',catFile,
	                  '-c',os.path.join(configDir,'boksolve.scamp')]
	def add_scamp_pars(scamp_pars):
		scamp_cmd = copy(scamp_cmd_base)
		for k,v in scamp_pars.items():
			scamp_cmd.extend(['-'+k,str(v)])
		return scamp_cmd
	try:
		os.unlink(wcsFile)
	except:
		pass
	#
	# FIRST PASS
	#
	scamp_pars = {
	  'DISTORT_DEGREES':1,
	  'ASTRCLIP_NSIGMA':3.0,
	  'CHECKPLOT_TYPE':'NONE',
	}
	if verbose > 0:
		scamp_pars['VERBOSE_TYPE'] = 'FULL'
	if refStarCatFile is not None and os.path.exists(refStarCatFile):
		scamp_pars['ASTREFCAT_NAME'] = refStarCatFile
		scamp_pars['ASTREF_CATALOG'] = 'FILE'
	else:
		scamp_pars['ASTREF_CATALOG'] = 'SDSS-R9'
		scamp_pars['ASTREF_BAND'] = filt
		scamp_pars['SAVE_REFCATALOG'] = 'Y'
		# see below
		#scamp_pars['ASTREFCAT_NAME'] = os.path.basename(refStarCatFile)
		#scamp_pars['REFOUT_CATPATH'] = refCatPath
	scamp_cmd = add_scamp_pars(scamp_pars)
	if verbose > 1:
		print ' '.join(scamp_cmd)
	elif verbose > 0:
		print 'first pass scamp_solve for ',imageFile
	rv = subprocess.call(scamp_cmd)
	tmpAhead = catFile.replace('.fits','.ahead') 
	shutil.move(headf,tmpAhead)
	if refStarCatFile is not None and scamp_pars['ASTREF_CATALOG'] != 'FILE':
		# scamp automatically names the cached reference file, and as far
		# as I can tell ignores the value of ASTREFCAT_NAME
		# take the auto-saved file and rename it
		tmpfn = min(glob.iglob('%s_*.cat'%scamp_pars['ASTREF_CATALOG']),
		            key=os.path.getctime)
		if verbose > 0:
			print 'tmpfn=',tmpfn
		#shutil.move(tmpfn,refStarCatFile)
		# XXX for some reason RA/Dec are being flipped
		f = fits.open(tmpfn)
		dec = f[2].data['X_WORLD'].copy()
		f[2].data['X_WORLD'][:] = f[2].data['Y_WORLD']
		f[2].data['Y_WORLD'][:] = dec
		f.writeto(refStarCatFile)
		os.unlink(tmpfn)
	#
	# SECOND PASS
	#
	if refStarCatFile is not None:
		scamp_pars['ASTREF_CATALOG'] = 'FILE'
		scamp_pars['ASTREFCAT_NAME'] = refStarCatFile
	try:
		del scamp_pars['ASTREF_BAND']
	except:
		pass
	scamp_pars['SAVE_REFCATALOG'] = 'N'
	scamp_pars['POSITION_MAXERR'] = 0.25
	scamp_pars['POSANGLE_MAXERR'] = 0.5
	scamp_pars['CROSSID_RADIUS'] = 5.0
	scamp_pars['DISTORT_DEGREES'] = 3
	scamp_pars['MOSAIC_TYPE'] = 'FIX_FOCALPLANE'
	if check_plots:
		del scamp_pars['CHECKPLOT_TYPE']
	scamp_cmd = add_scamp_pars(scamp_pars)
	if verbose > 1:
		print ' '.join(scamp_cmd)
	elif verbose > 0:
		print 'second pass scamp_solve for ',imageFile
	rv = subprocess.call(scamp_cmd)
	shutil.move(headf,wcsFile)
	os.unlink(tmpAhead)
	#
	if savewcs:
		configFile = os.path.join(configDir,'wcsput.missfits')
		missfits_cmd = ['missfits','-c',configFile,imageFile]
		if verbose > 1:
			print ' '.join(missfits_cmd)
		rv = subprocess.call(missfits_cmd)

