#!/usr/bin/env python

import os,shutil
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits

from bokpipe.badpixels import build_mask_from_flat
from bokpipe import bokpl,bokobsdb
from bokpipe import __version__ as pipeVersion

def set_bass_defaults(args):
	if args.rawdir is None:
		args.rawdir = os.environ['BASSDATA']
	if args.output is None:
		outdir = os.path.join(os.environ['BASSRDXDIR'],'reduced',pipeVersion)
		args.output = outdir
	if args.obsdb is None:
		args.obsdb = os.path.join('config','bassobsdb.fits')
	return args

def config_bass_data(dataMap,args,season=None):
	if args.band is None:
		dataMap.setFilters(['g','bokr'])
	dataMap.setFringeFilters(['bokr'])
	if season:
		sfx = '_%s' % season
	else:
		sfx = ''
	if args.badpixdir:
		for sfx2 in ['','_x4']:
			bpfn = 'BadPixMask%s%s.fits' % (sfx,sfx2)
			bpfile = os.path.join(args.badpixdir,bpfn)
			outfile = os.path.join(dataMap.getCalDir(),bpfn)
			if not os.path.exists(outfile):
				print 'copying ',bpfile,' to ',outfile
				shutil.copy(bpfile,outfile)
	dataMap.setCalMap('badpix','master',fileName='BadPixMask%s'%sfx)
	dataMap.setCalMap('badpix4','master',fileName='BadPixMask%s_x4'%sfx)
	dataMap.setCalMap('ramp','master',fileName='BiasRamp%s'%sfx)
	return dataMap

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

def get_observing_season(dataMap):
	year = np.unique([ utd[:4] for utd in dataMap.getUtDates() ])
	if len(year) > 1:
		raise ValueError("Currently processing only one observing season "
		                 "at a time is supported")
	return year[0]

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
		elif keep.sum() > self.maxNImg:
			jj = np.where(keep)[0]
			np.random.shuffle(jj)
			keep[jj[self.maxNImg:]] = False
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
	parser.add_argument('--makebpmask',type=str,
	                help='make quick badpix mask from flat <FILENAME>')
	parser.add_argument('--badpixdir',type=str,
	                help='location of bad pixel mask masters')
	args = parser.parse_args()
	args = set_bass_defaults(args)
	dataMap = bokpl.init_data_map(args)
	season = get_observing_season(dataMap)
	dataMap = config_bass_data(dataMap,args,season)
	kwargs = {}
	kwargs['illum_filter_fun'] = IllumFilter()
	if args.finish:
		finish_up(args,dataMap)
	elif args.makebpmask:
		build_mask_from_flat(dataMap('cal')(args.makebpmask),
		                     dataMap.getCalMap('badpix').getFileName(),
		                     dataMap.getCalDir(),
		                     verbose=args.verbose,
		                     debug=args.debug)
	else:
		bokpl.run_pipe(dataMap,args,**kwargs)

