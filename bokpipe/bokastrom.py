#!/usr/bin/env python

import os
import glob
import subprocess
from copy import copy
import numpy as np
import fitsio

def scamp_solve(imageFile,catFile,wcsFile,refStarCatFile=None,
                filt='g',savewcs=False,overwrite=False):
	configDir = os.path.join(os.path.split(__file__)[0],'config')
	headf = catFile.replace('.fits','.head') # named automatically by scamp?
	if not wcsFile.endswith('.ahead'):
		wcsFile += '.ahead'
	if not overwrite and os.path.exists(wcsFile):
		print wcsFile,' already exists, skipping'
		return
	#
#	rmfield = logs[utdate]['objectName'][frame]
#	filt = logs[utdate]['filter'][frame]
#	astrefcatfn = '%s_ccd%d_%s_ref.cat' % (rmfield,ccdNum,filt)
#	astrefcatf = os.path.join(wcscatdir,astrefcatfn)
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
	if refStarCatFile is not None and os.path.exists(refStarCatFile):
		scamp_pars['ASTREFCAT_NAME'] = os.path.basename(refStarCatFile)
		scamp_pars['ASTREF_CATALOG'] = 'FILE'
		scamp_pars['REFOUT_CATPATH'] = os.path.dirname(refStarCatFile)
	else:
		scamp_pars['ASTREF_CATALOG'] = 'SDSS-R9'
		scamp_pars['ASTREF_BAND'] = filt
		scamp_pars['SAVE_REFCATALOG'] = 'Y'
	scamp_cmd = add_scamp_pars(scamp_pars)
	print ' '.join(scamp_cmd)
	rv = subprocess.call(scamp_cmd)
	os.rename(headf,wcsFile)
	if False: #scamp_pars['ASTREF_CATALOG'] != 'FILE':
		# scamp automatically names the cached reference file, and as far
		# as I can tell ignores the value of ASTREFCAT_NAME
		# take the auto-saved file and rename it
		tmpfn = min(glob.iglob('%s_*.cat'%scamp_pars['ASTREF_CATALOG']),
		            key=os.path.getctime)
		print 'tmpfn=',tmpfn
		os.rename(tmpfn,astrefcatf)
	#
	# SECOND PASS
	#
	scamp_pars['ASTREF_CATALOG'] = 'FILE'
	scamp_pars['REFOUT_CATPATH'] = os.path.dirname(refStarCatFile)
	try:
		del scamp_pars['ASTREF_BAND']
	except:
		pass
	scamp_pars['SAVE_REFCATALOG'] = 'N'
	scamp_pars['POSITION_MAXERR'] = 0.1
	scamp_pars['POSANGLE_MAXERR'] = 0.5
	scamp_pars['CROSSID_RADIUS'] = 2.0
	scamp_pars['DISTORT_DEGREES'] = 3
	scamp_cmd = add_scamp_pars(scamp_pars)
	print ' '.join(scamp_cmd)
	rv = subprocess.call(scamp_cmd)
	os.rename(headf,wcsFile)
	#
	if savewcs:
		# XXX move this to config
		missfits_cmd = ['missfits','-c',
	       os.path.join(os.environ['BOK90PRIMEDIR'],'pro','default.missfits'),
		                imageFile]
		rv = subprocess.call(missfits_cmd)

