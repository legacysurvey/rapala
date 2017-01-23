#!/usr/bin/env python

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import argparse

from bokpipe.bokmkimage import make_fov_image_fromfile

parser = argparse.ArgumentParser()
parser.add_argument("fitsFile",type=str,
                    help="input FITS image")
parser.add_argument("imgFile",type=str,nargs='?',
                    help="output image file")
parser.add_argument("--mask",type=str,
                    help="input mask image")
parser.add_argument("--nbin",type=int,default=1,
                    help="output image file")
parser.add_argument("--coordsys",type=str,default='sky',
                    help="coordinate system")
parser.add_argument("--stretch",type=str,default='linear',
                    help="stretch ([linear])")
parser.add_argument("--interval",type=str,default='zscale',
                    help="interval\n"+
                         " zscale (default)\n"+
                         " rms (combine with range)\n"+
                         " fixed (combine with range)")
parser.add_argument("--imrange",type=str,
                    help="range (vmin,vmax), (siglo,sighi)")
parser.add_argument("--contrast",type=float,default=0.25,
                    help="contrast for zscale [default: 0.25]")
parser.add_argument("--cmap",type=str,
                    help="color map")
args = parser.parse_args()

if args.imgFile:
	outFile = args.imgFile
else:
	outFile = args.fitsFile.replace('.fits','.png')

try:
	imrange = [float(v) for v in args.imrange.split(',')]
except:
	imrange = None

make_fov_image_fromfile(args.fitsFile,outFile,mask=args.mask,
                        nbin=args.nbin,coordsys=args.coordsys,
                        stretch=args.stretch,interval=args.interval,
                        imrange=imrange,contrast=args.contrast,
                        cmap=args.cmap)

