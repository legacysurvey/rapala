#!/usr/bin/env python

import os
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits

from bokpipe import bokpl,bokobsdb
from bokpipe import __version__ as pipeVersion

def set_bass_defaults(args):
	if args.rawdir is None:
		args.rawdir = os.environ['BASSDATA']
	if args.output is None:
		outdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',pipeVersion)
		args.output = outdir
	if args.obsdb is None:
		args.obsdb = os.path.join('config','basslog.fits')
	return args

def finish_up(args,dataMap):
	'''rename files and move them to a separate directory'''
	import shutil
	from bokpipe.bokastrom import put_wcs
	indir = args.output
	outdir = os.path.join(indir,'nov15data')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	for filt in dataMap.iterFilters():
		filesAndNames = dataMap.getFiles(with_objnames=True)
		outdirb = os.path.join(outdir,filt[-1])
		if not os.path.exists(outdirb):
			os.mkdir(outdirb)
		for f,n in zip(*filesAndNames):
			n = n.replace('bokr','r')
			outf = os.path.join(outdirb,n)
			print f,outf
			for ext in ['.fits','.wht.fits','.psf.fits','.ahead','.cat.fits']:
				shutil.copy(os.path.join(indir,f)+ext,outf+ext)
			put_wcs(outf+'.fits')
			os.remove(outf+'.ahead')

class IllumFilter(object):
	minNImg = 15
	maxNImg = 30
	maxCounts = 10000
	def __call__(self,obsDb,ii):
		keep = np.ones(len(ii),dtype=bool)
		# this whole night is bad due to passing clouds
		#keep[:] ^= obsDb['utDate'][ii] == '20140612'
		if keep.sum()==0:
			return keep
		elif keep.sum() > maxNImg:
			jj = np.where(keep)[0]
			np.random.shuffle(jj)
			keep[jj[self.maxNimg:]] = False
		print 'selected %d images for illumination correction' % keep.sum()
		return keep

if __name__=='__main__':
	import sys
	import argparse
	parser = argparse.ArgumentParser()
	parser = bokpl.init_file_args(parser)
	parser = bokpl.init_pipeline_args(parser)
	parser.add_argument('--finish',action='store_true',
	                help='rename files')
	args = parser.parse_args()
	args = set_bass_defaults(args)
	dataMap = bokpl.init_data_map(args)
	if args.band is None:
		dataMap.setFilters(['g','bokr'])
	dataMap = bokpl.set_master_cals(dataMap)
	kwargs = {}
	kwargs['illum_filter_fun'] = IllumFilter()
	if args.finish:
		finish_up(args,dataMap)
	else:
		bokpl.run_pipe(dataMap,args,**kwargs)

