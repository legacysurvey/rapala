#!/usr/bin/env python

import sys

from bokpipe import *

if __name__=='__main__':
	outputFile = sys.argv[1]
	inputFiles = sys.argv[2:]
	method = 'mean'
	stackPars = {}
	stackPars['scale'] = 'normalize_median'
	stackPars['reject'] = 'sigma_clip'
	stackPars['nsplit'] = 8
	if method == 'mean':
		stackFun = bokutil.ClippedMeanStack(**stackPars)
	elif method == 'median':
		stackFun = bokutil.MedianStack(**stackPars)
	_inputFiles = []
	# accept an IRAF-style list from a file if args starts with '@'
	for f in inputFiles:
		if f.startswith('@'):
			with open(f[1:]) as _f:
				_inputFiles.extend(_f.read().splitlines())
		else:
			_inputFiles.append(f)
	print 'input: ',_inputFiles
	print 'output: ',outputFile
	stackFun.stack(inputFiles,outputFile)

