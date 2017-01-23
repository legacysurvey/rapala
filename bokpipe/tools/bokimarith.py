#!/usr/bin/env python

import sys
from bokpipe import bokproc

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input",nargs='+',type=str,
                    help="input FITS image(s)")
parser.add_argument("operator",type=str,
                    help="mathematical operator (+-x/)")
parser.add_argument("operand",type=str,
                    help="image or value used as operand")
parser.add_argument("-o","--output",type=str,
                    help="output image name (default is same)")
parser.add_argument("-p","--prefix",type=str,default='',
                    help="prefix for output file names")
parser.add_argument("-s","--suffix",type=str,default='',
                    help="suffix for output file names")
args = parser.parse_args()

def renamemap(f):
	return args.prefix+f.replace('.fits','')+args.suffix+'.fits'

if args.output:
	namemap = lambda f: args.output
else:
	namemap = renamemap

imarith = bokproc.BokImArith(args.operator,args.operand,output_map=namemap)

imarith.process_files(args.input)


