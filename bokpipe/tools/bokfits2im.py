#!/usr/bin/env python

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
parser.add_argument("--vmin",type=float,
                    help="minimum range")
parser.add_argument("--vmax",type=float,
                    help="maximum range")
parser.add_argument("--siglo",type=float,default=2.5,
                    help="minimum range")
parser.add_argument("--sighi",type=float,default=5.0,
                    help="maximum range")
parser.add_argument("--cmap",type=str,
                    help="color map")
args = parser.parse_args()

if args.imgFile:
	outFile = args.imgFile
else:
	outFile = args.fitsFile.replace('.fits','.png')

make_fov_image_fromfile(args.fitsFile,outFile,mask=args.mask,
                        nbin=args.nbin,coordsys=args.coordsys,
                        vmin=args.vmin,vmax=args.vmax,
                        lo=args.siglo,hi=args.sighi,
                        cmap=args.cmap)

