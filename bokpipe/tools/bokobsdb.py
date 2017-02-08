#!/usr/bin/env python

import numpy as np
from astropy.table import Table

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("obsdb",type=str,
                    help="input observations database")
parser.add_argument('-b','--band',type=str,default=None,
                help='band to include [default=both]')
parser.add_argument('--frames',type=str,default=None,
                help='frames to include (i1,i2) [default=all]')
parser.add_argument('-t','--imtype',type=str,default=None,
                help='specify image type to include')
parser.add_argument('-u','--utdate',type=str,default=None,
                help='UT date(s) to include [default=all]')
parser.add_argument('--night',type=str,default=None,
                help='night(s) to include [default=all]')
parser.add_argument('--fields',type=str,
                help='additional fields to display')
args = parser.parse_args()

fields = ['frameIndex','fileName','imType',
          'filter','objName','expTime']

if args.fields:
	fields += args.fields.split(',')

obsdb = Table.read(args.obsdb)

frames = np.ones(len(obsdb),dtype=bool)
if args.band:
	frames &= np.in1d(obsdb['filter'],args.band)
if args.imtype:
	frames &= np.in1d(obsdb['imType'],args.imtype)
if args.utdate:
	frames &= np.in1d(obsdb['utDate'],args.utdate)

obsdb[frames][fields].more()

