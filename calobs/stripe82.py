#!/usr/bin/env python

import sys
from astropy import units as u
from astropy.coordinates import SkyCoord

ra = float(sys.argv[1])

filt = 'g'
#filt = 'bokr'
exptime = 100.

dec = 0.0
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

outf = open("s82cal_%s_ra%d.txt"%(filt,ra),'w')

outf.write('movefilter %s\n' % filt)

for i,(ra_off,dec_off) in enumerate(offsets):
	name = 's82cal%s_ra%d_%d' % (filt,ra,i+1)
	c = SkyCoord(ra=(ra+ra_off)*u.degree,dec=(dec+dec_off)*u.degree)
	s = c.to_string('hmsdms',sep='',precision=2)
	l = "obs %.1f object '%s' %d %s %s 2000.0" % \
	      (exptime,name,nexp,filt,s)
	print l
	outf.write(l+"\n")

outf.close()

