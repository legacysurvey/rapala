#!/usr/bin/env python

import os,sys
import numpy as np
import fitsio
from bokpipe.bokutil import BokMefImage,array_stats,stats_region
from bokpipe import bokproc

import argparse

class CCDMedianBackgroundFit(bokproc.BackgroundFit):
	def __init__(self,fits,**kwargs):
		super(CCDMedianBackgroundFit,self).__init__(fits,**kwargs)
		self.statsReg = stats_region('ccd_central_quadrant')
		self.medBack = {}
		for extn,data,hdr in fits:
			self.medBack[extn] = array_stats(data[self.statsReg],
			                                 method='median',clip=True)
			self.imShape = data.shape
	def get(self,extn):
		return self.medBack[extn] + np.zeros(self.imShape,dtype=np.float32)

parser = argparse.ArgumentParser()
parser.add_argument("input",type=str,
                    help="input FITS image")
parser.add_argument("-f","--fit",type=str,default='spline',
                    help="fit method (poly|[spline])")
parser.add_argument("--order",type=int,default=1,
                    help="fit order [default=1]")
parser.add_argument("--knots",type=int,default=2,
                    help="number of spline knots [default=2]")
parser.add_argument("--nbin",type=int,default=64,
                    help="amount of binning to apply before fit[default=64]")
parser.add_argument("-o","--output",type=str,
                    help="output FITS image")
parser.add_argument("-m","--mask",type=str,
                    help="mask FITS image")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
args = parser.parse_args()

fileName = args.input

fits = BokMefImage(fileName,mask_file=args.mask,read_only=True)

if args.fit == 'spline':
	backfit = bokproc.SplineBackgroundFit(fits,nKnots=args.knots,
	                                      order=args.order,nbin=args.nbin)
elif args.fit == 'poly':
	backfit = bokproc.PolynomialBackgroundFit(fits,
	                                          order=args.order,nbin=args.nbin)
elif args.fit == 'ccdmedian':
	backfit = CCDMedianBackgroundFit(fits)
else:
	raise ValueError

if not args.output:
	outfn = fileName.replace('.fits','_back.fits')
else:
	outfn = args.output

outFits = fitsio.FITS(outfn,'rw',clobber=True)

hdr = fits.get_header(0)
outFits.write(None,header=hdr)

for extn,data,hdr in fits:
	outFits.write(backfit.get(extn),extname=extn,header=hdr)

outFits.close()


