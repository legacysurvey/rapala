#!/usr/bin/env python

import argparse

from bokpipe import bokobsdb

parser = argparse.ArgumentParser()
parser.add_argument("inputDirs",type=str,nargs='+',
                    help="directories containing FITS images")
parser.add_argument("-o","--output",type=str,default="boklog.fits",
                    help="output log file")
parser.add_argument("-f","--filters",type=str,
                    help="filter list")
args = parser.parse_args()

bokobsdb.generate_log(args.inputDirs,args.output,
                      filters=args.filters,
                      objFilter=None,
                      filePattern=None)

