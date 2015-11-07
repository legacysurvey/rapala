#!/usr/bin/env python

import os,sys

from bokpipe.bokproc import combine_ccds

if __name__=='__main__':
	inputFile = sys.argv[1]
	if len(sys.argv)>2:
		outputFile = sys.argv[2]
	else:
		outputFile = inputFile.replace('.fits','_ccd4.fits')
	kwargs = {}
	combine_ccds([inputFile,],output_map=lambda s: outputFile,**kwargs)

