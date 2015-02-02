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
	elevation = 90 - np.degrees(np.arccos(1/airmass))
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

