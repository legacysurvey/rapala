#!/usr/bin/env/python

import argparse

from bokpipe.bokextract import sextract

parser = argparse.ArgumentParser()
parser.add_argument("image",type=str,
                    help="input FITS image")
parser.add_argument("catalog",type=str,
                    help="output FITS catalog")
parser.add_argument("-r","--redo",action="store_true",
                    help="redo (overwrite existing")
parser.add_argument('-v','--verbose',action='count',
                    help='increase output verbosity')
args = parser.parse_args()

sextract(args.image,args.catalog,overwrite=args.redo,verbose=args.verbose)

