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
	print 'input: ',inputFiles
	print 'output: ',outputFile
	stackFun.stack(inputFiles,outputFile)

