#!/usr/bin/env python

import os
import argparse
import fitsio

from bokpipe.bokproc import combine_ccds

parser = argparse.ArgumentParser()
parser.add_argument("input",type=str,
                    help="input FITS image")
parser.add_argument("-o","--output",type=str,
                    help="output FITS image")
parser.add_argument("--split",action="store_true",
                    help="split into per-ccd files (instead of join)")
parser.add_argument("--flatnorm",action="store_true",
                    help="normalize to unity at CCD centers for flats")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
args = parser.parse_args()

if args.split:

	f = fitsio.FITS(args.input)
	for ccd,hdu in enumerate(f[1:],start=1):
		if args.output is None:
			outputFile = args.input.replace('.fits','.ccd%d.fits'%ccd)
		else:
			outputFile = args.output.replace('.fits','.ccd%d.fits'%ccd)
		outputFits = fitsio.FITS(outputFile,'rw')
		outputFits.write(hdu.read(),header=hdu.read_header())
		outputFits.close()

else:

	if args.output is None:
		outputFile = args.input.replace('.fits','_4ccd.fits')
	else:
		outputFile = args.output
	
	kwargs = {'clobber':True,'apply_flat_norm':args.flatnorm}
	
	combine_ccds([args.input,],output_map=lambda s: outputFile,**kwargs)

