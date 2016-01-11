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

if __name__=='__main__':
	import sys
	import argparse
	parser = argparse.ArgumentParser()
	parser = bokpl.init_file_args(parser)
	parser = bokpl.init_pipeline_args(parser)
	args = parser.parse_args()
	args = set_bass_defaults(args)
	dataMap = bokpl.init_data_map(args)
	dataMap = bokpl.set_master_cals(dataMap)
	bokpl.run_pipe(dataMap,args)

