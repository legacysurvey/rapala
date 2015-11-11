#!/usr/bin/env python

import sys
from astropy import units as u
from astropy.coordinates import SkyCoord

dec = 0.0

try:
	name = sys.argv[1]
	ra = float(sys.argv[2])
	filt = sys.argv[3]
	assert filt in ['g','bokr']
	exptime = float(sys.argv[4])
	if len(sys.argv) > 5:
		dec = float(sys.argv[5])
except:
	print 'Usage: caldither.py name ' \
	      'ra(deg) filter(g|bokr) exptime(s) [dec(deg)=0]'
	sys.exit(1)

nexp = 1

ra_gap = 170./3600
dec_gap = 55./3600

offsets = [
  (0,0),
  (0.25,0),
  (0.25,0.25),
  (0,0.25),
  (0,-0.5-dec_gap),
  (0.5+ra_gap,0.0),
]

outf = open("%s_%s.txt"%(name,filt),'w')

outf.write('movefilter %s\n' % filt)

for i,(ra_off,dec_off) in enumerate(offsets):
	pname = '%s%s_%d' % (name,filt,i+1)
	c = SkyCoord(ra=(ra+ra_off)*u.degree,dec=(dec+dec_off)*u.degree)
	s = c.to_string('hmsdms',sep='',precision=2)
	l = "obs %.1f object '%s' %d %s %s 2000.0" % \
	      (exptime,pname,nexp,filt,s)
	print l
	outf.write(l+"\n")

outf.close()

