#!/usr/bin/env python

import os
import glob
import numpy as np
from astropy.table import Table

from bokpipe import bokpl,bokobsdb
from bokpipe import __version__ as pipeVersion

def set_bass_defaults(args):
	if args.rawdir is None:
		args.rawdir = os.environ['BASSDATA']
	if args.output is None:
		args.output = os.path.join('reduced',pipeVersion)
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
		outdirb = os.path.join(outdir,filt)
		if not os.path.exists(outdirb):
			os.mkdir(outdirb)
		for f,n in zip(*filesAndNames):
			outf = os.path.join(outdirb,n)
			print f,outf
			for ext in ['.fits','.psf.fits','.ahead','.cat.fits']:
				shutil.copy(os.path.join(indir,f)+ext,outf+ext)
			put_wcs(outf+'.fits')
			os.remove(outf+'.ahead')

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
	dataMap = bokpl.set_master_cals(dataMap)
	if args.finish:
		finish_up(args,dataMap)
	else:
		bokpl.run_pipe(dataMap,args)

