#!/usr/bin/env python

import os,sys
import numpy as np
import fitsio
from bokpipe.bokutil import BokMefImage,array_stats,stats_region
from bokpipe import bokproc

import argparse

class CCDMedianBackgroundFit(bokproc.BackgroundFit):
	def __init__(self,fits,**kwargs):
		self.statsReg = kwargs.pop('stats_reg','ccd_central_quadrant')
		self.method = kwargs.pop('method','median')
		super(CCDMedianBackgroundFit,self).__init__(fits,**kwargs)
		self.statsPix = stats_region(self.statsReg)
		self.medBack = {}
		for extn,data,hdr in fits:
			self.medBack[extn] = array_stats(data[self.statsPix],
			                                 method=self.method,clip=True)
			self.imShape = data.shape
	def get(self,extn):
		return self.medBack[extn] + np.zeros(self.imShape,dtype=np.float32)

parser = argparse.ArgumentParser()
parser.add_argument("input",type=str,
                    help="input FITS image")
parser.add_argument("-f","--fit",type=str,default='spline',
                    help="fit method (poly|[spline]|ccdmedian|ccdmean|ccdmode)")
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
parser.add_argument("-s","--statsreg",type=str,default='ccd_central_quadrant',
                    help="region to use for statistics (mean/median/mode)")
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
	backfit = CCDMedianBackgroundFit(fits,stats_reg=args.statsreg,
	                                 method=args.fit[3:])
else:
	raise ValueError

if not args.output:
	outfn = fileName.replace('.fits','_back.fits')
else:
	outfn = args.output

backfit.write(outfn,clobber=True)

