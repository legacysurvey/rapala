#!/usr/bin/env python

import os
import re
import glob
from copy import copy
import multiprocessing
import numpy as np
from numpy.core.defchararray import add as char_add
import fitsio
from astropy.table import Table

from bokpipe import bokpl
from bokpipe import __version__ as pipeVersion

def load_darksky_frames(filt):
	darkSkyFrames = np.loadtxt(os.path.join('config',
	                                        'bokrm_darksky_%s.txt'%filt),
	                           dtype=[('utDate','S8'),('fileName','S35'),
	                                  ('skyVal','f4')])
	# a quick pruning of the repeat images
	return darkSkyFrames[::2]

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser = bokpl.init_file_args(parser)
	parser = bokpl.init_pipeline_args(parser)
	parser.add_argument('--darkskyframes',action='store_true',
	                help='load only the dark sky frames')
	args = parser.parse_args()
	if args.rawdir is None:
		args.rawdir = os.environ['BOK90PRIMERAWDIR']
	if args.output is None:
		args.output = os.path.join(os.environ['BOK90PRIMEOUTDIR'],
		                           pipeVersion)
	if args.obsdb is None:
		args.obsdb = os.path.join('config','sdssrm-bok2014.fits')
	dataMap = bokpl.init_data_map(args)
	dataMap = bokpl.set_master_cals(dataMap)
	if args.darkskyframes:
		# must select a band
		if bands is None or bands not in ['g','i']:
			raise ValueError("Must select a band for dark sky frames (-b)")
		frames = load_darksky_frames(bands)
		dataMap.setFrameList(frames['utDate'],frames['fileName'])
	bokpl.run_pipe(dataMap,args)

