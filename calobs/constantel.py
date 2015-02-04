#!/usr/bin/env python

import itertools
import numpy as np
from scipy.interpolate import interp1d
import pidly

import matplotlib.pyplot as plt

def calc_az_track(utdate,utStart,maxDuration,exptime,overhead,elevation,
                  offsetScale=3.0,jumpFrac=0.1,
                  idl=None):
	'''calls the idl routine for generating a fixed-elevation track.
	     utdate: string, like '20150101'
	     utStart: time in UT hours to start the track (e.g., 2.0 == 02:00UT)
	     maxDuration: maximum duration of track in hours
	     exptime: individual exposure time in seconds
	     overhead: overhead for individual exposure in seconds
	     elevation: the elevation of the track in degrees
	     offsetScale: the scale of random AZ offsets in degrees
	     jumpFrac: fraction of offsets to be in negative AZ
	'''
	# the time between exposures, including overheads
	shotTime = (exptime+overhead)/3600.
	if idl is None:
		idl = pidly.IDL()
	idl.year = int(utdate[:4])
	idl.month = int(utdate[4:6])
	idl.day = int(utdate[6:8])
	idl.hours = np.arange(utStart,utStart+maxDuration,shotTime)
	idl.elevation = elevation
	#seed = 1 #np.random.randint(2**63-1)
	idl("aztrack,year,month,day,hours,elevation,ra,dec,az,utout,"+
	    "jumpFrac=%f,offsetScale=%f" % (jumpFrac,offsetScale))
	    #"jumpFrac=%f,offsetScale=%f,seed=%d" % (jumpFrac,offsetScale,seed))
	if np.all(idl.ra < 0):
		return dict(duration=0)
	duration = (idl.utout[-1]+shotTime) - idl.utout[0]
	return dict(ra=idl.ra,dec=idl.dec,ut=idl.utout,az=idl.az,
	            duration=duration)

def airmass2el(airmass):
	return 90 - np.degrees(np.arccos(1/airmass))

def plot_season(airmass=1.4,**kwargs):
	'''Plot the pointings for the bright-time calibration strategy at
	   fixed airmass using bright time observing dates in 2015.
	   kwargs:
	     airmass: fixed airmass for tracks (default 1.4)
	     minDuration: which is the minimum time for a track in hours 
	                  (i.e., don't include tracks shorter than minDuration)
	                  the default is 0.25 == 15 minutes
	     utStart: UT time at start of track (default is 2.0==02:00UT==19:00MST)
	     utEnd:   UT time at end of track (default is 13.0==13:00UT==06:00MST)
	     utBreak: UT times to "take a break" (e.g., for other observing)
	               default is (5.0,10.0) which assumes 5 hours in the middle
	               of the night will be used for something else
	               set to (np.inf,None) to never take a break
	     additional kwargs go to calc_az_track()
	'''
	utdates = ['20150203','20150305','20150402','20150505','20150706']
	elevation = airmass2el(airmass)
	minDuration = kwargs.get('minDuration',0.25) # 15 minutes
	utStart = kwargs.get('utStart',2.0)
	utEnd = kwargs.get('utEnd',13.0)
	utBreak = kwargs.get('utBreak',(5.0,10.0))
	alltracks = {}
	idl = pidly.IDL()
	for utdate,c in zip(utdates,itertools.cycle('rgbcymk')):
		alltracks[utdate] = []
		uttime = utStart
		utskip = True
		while uttime < utEnd:
			# random duration between 1/2 hour and 1.5 hours
			duration = 0.5 + np.random.rand()
			track = calc_az_track(utdate,uttime,duration,50.,50.,elevation,
			                      idl=idl,**kwargs)
			if track['duration'] < minDuration:
				# too short, try again
				continue
			print '%s %7.2f %7.2f' % (utdate,uttime,track['duration']*60)
			plt.plot(track['ra'],track['dec'],'%ss-'%c)
			uttime += track['duration']
			alltracks[utdate].append(track)
			if uttime > utBreak[0] and utskip:
				# skip a chunk in the middle of the night
				uttime = utBreak[1]
				utskip = False
	for ut in utdates:
		for j,t in enumerate(alltracks[ut]):
			print '%s %3d %7.2f %7.2f %7.2f' % \
			       (ut,j,t['dec'][0],t['dec'][-1],t['duration']*60)
	return alltracks

def formatradec(_ra,_dec,withdelim=False):
	from astrolib import coords
	p = coords.Position((_ra,_dec))
	if withdelim:
		return p.hmsdms()
	s_ra,s_dec = p.hmsdms().replace(':','').split()
	# remove one/two digit of precision
	return (s_ra[:-1],s_dec[:-2])

def formatut(ut,full=False):
	uth = int(ut)
	utm = int(60 * (ut % uth))
	if full:
		if utm==0:
			uts = 0
		else:
			uts = int(60 * (60*(ut-uth) % utm))
		return '%02d:%02d:%02d' % (uth,utm,uts)
	else:
		return '%02d%02d' % (uth,utm)

def dump_track(track,utdate,exposureTime,filt):
	from itertools import count
	utstr = formatut(track['ut'][0])
	scriptf = open('basscal_%s_%s_%s.txt' % (utdate,utstr,filt),'w')
	logf = open('basscal_%s_%s_%s.log' % (utdate,utstr,filt),'w')
	for i,ut,ra,dec in zip(count(),track['ut'],track['ra'],track['dec']):
		ras,decs = formatradec(ra,dec)
		imtitle = 'cal%s%s%s_%02d' % (utdate,utstr,filt,i)
		scriptf.write("obs %.1f object '%s' 1 %s %s %s 2000.0\r\n" % 
		              (exposureTime,imtitle,filt,ras,decs))
		coo = formatradec(ra,dec,True)
		logf.write('%s   %s   %s\n' % (formatut(ut,True),imtitle,coo))
	scriptf.close()
	logf.close()

def make_track(utdate,uttime,**kwargs):
	airmass = kwargs.get('airmass',1.4)
	exposureTime = kwargs.get('exposureTime',50.)
	overheadTime = kwargs.get('overheadTime',50.)
	duration = kwargs.get('duration',1.0)
	offsetScale = kwargs.get('offsetScale',3.0)
	jumpFrac = kwargs.get('jumpFrac',0.1)
	filt = kwargs.get('filt','g')
	elevation = airmass2el(airmass)
	uttime = float(uttime)
	print
	print 'Constructing track at airmass=%.2f' % airmass
	print '  Duration: %.1f hr with exposure time = %.1fs (overhead %.1fs)' % \
	            (duration,exposureTime,overheadTime)
	print '  random params: offsetScale=%.1f deg, jumpFrac=%.1f' % \
	            (offsetScale,jumpFrac)
	print '  elevation is %.1f degrees' % elevation
	print
	track = calc_az_track(utdate,uttime,duration,exposureTime,overheadTime,
	                      elevation,offsetScale=offsetScale,jumpFrac=jumpFrac)
	if 'ut' not in track:
		return None
	dump_track(track,utdate,exposureTime,filt)
	return track

if __name__=='__main__':
	import sys,getopt
	try:
		opts,args = getopt.getopt(sys.argv[1:],
		                          "a:d:e:o:",
		                         ['airmass=','duration=',
		                          'exposureTime=','overheadTime=',])
	except getopt.GetoptError as err:
		print str(err)
		sys.exit(2)
	if len(args) != 2:
		print 'must provide utdate and time as arguments'
		sys.exit(2)
	kwargs = {}
	for k,v in opts:
		kwargs[k.lstrip('--')] = float(v)
	ntry = 0
	while ntry < 10:
		t = make_track(*args,**kwargs)
		if t is not None:
			break
		ntry += 1

