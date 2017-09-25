#!/usr/bin/env python

import os,sys
import time
from collections import defaultdict
import numpy as np
from astropy.table import Table
from bokpipe.bokdm import load_caldb,write_caldb

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("caldb",type=str,
                    help="input calibrations database")
parser.add_argument("--obsdb",type=str,
                    help="input observations database")
parser.add_argument('-b','--band',type=str,default=None,
                help='band to include [default=both]')
parser.add_argument('-t','--caltype',type=str,default=None,
                help='specify calibraion type [zero|flat|illum|skyflat|fringe]')
parser.add_argument('-u','--utdate',type=str,default=None,
                help='UT date(s) to include [default=all]')
parser.add_argument('--night',type=str,default=None,
                help='night(s) to include [default=all]')
parser.add_argument('-l','--listframes',action="store_true",
                help='list individual frames')
parser.add_argument('-d','--delete',action="store_true",
                help='delete matching entries')
parser.add_argument('-f','--force',action="store_true",
                help='force deletion (no prompt)')
args = parser.parse_args()

caldb = load_caldb(args.caldb)

if args.listframes:
	if not args.obsdb:
		print "Must provide observations database [--obsdb]"
		sys.exit(1)
	obsdb = Table.read(args.obsdb)

if args.caltype:
	caltypes = [args.caltype]
else:
	caltypes = 'zero|flat|illum|skyflat|fringe'.split('|')

fields = ['frameIndex','utDate','fileName','objName','expTime']

if args.delete:
	newdb = defaultdict(list)

print
for t in caltypes:
	print '[[ %s ]]' % t
	if t not in caldb:
		print '<<None>>'
		continue
	for fn,utd,mjd,filt,frames in caldb[t]:
		if ( (args.utdate and not utd.startswith(args.utdate)) or
		     (args.band and not filt==args.band) ):
			if args.delete:
				newdb[t].append( (fn,utd,mjd,filt,frames) )
		else:
			print "{0} {1} {2}  {3}".format(utd,filt,fn,len(frames))
			if args.listframes:
				obsdb[fields][frames].pprint(-1)
			if args.delete and not args.force:
				resp = raw_input('delete [y/n/q]? ').lower()
				if resp.startswith('y'):
					pass
				elif resp.startswith('n'):
					newdb[t].append( (fn,utd,mjd,filt,frames) )
				elif resp.startswith('q'):
					sys.exit(0)
	print

if args.delete:
	for t in caldb.keys():
		if t not in newdb:
			newdb[t] = caldb[t]
	newdb = dict(newdb)
	t = os.path.getmtime(args.caldb)
	sfx = time.strftime('%Y%m%d%H%M%S',time.gmtime(t))
	bkpf = args.caldb.replace('.pkl','_%s.pkl'%sfx)
	write_caldb(bkpf,caldb)
	write_caldb(args.caldb,newdb)

