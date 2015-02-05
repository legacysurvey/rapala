#!/usr/bin/env python

import itertools
import numpy as np
from scipy.interpolate import interp1d
from astropy.time import Time
import astropy.coordinates as coo
from astropy import units as u

import matplotlib.pyplot as plt

# Generate a track of exposures at constant elevation. The track moves
# in increasing azimuth, except for a random frequency of negative
# azimuth moves.
#
# INPUT
#  year,month,day: date of observation
#  hours: array of UT times in hours for each exposure
#     e.g., hours=[1.0,1.5,2.0] means the exposures start at 
#              01:00UT,01:30UT,02:00UT
#  elevation: the constant elevation for the track
# OUTPUT
#  ra,dec,az,utout: vectors containing output tracks
#                   utout is a (possibly) truncated version of 'hours',
#                   in that the sequence may reach a limit in ra/dec 
# PARAMETERS
#  jumpFrac: fraction of exposures which get a negative AZ offset
#            default is 0.1 (10%). The negative offsets are at 3*offsetScale
#  offsetScale: maximum AZ offset between each exposure, in degrees
#            default is 3
#  decMin: The minimum declination allowed (default is 30), in degrees
#  decMax:     maximum   "             "   (default is 85)
#  raMin,raMax: same for RA, defaults are 30,300
#  gbMin: The minimum galactic latitude allowed (default is 17)
#  seed: fix the random seed for the track
#
def aztrack(startTime,shotTime,duration,fixedElevation,
            jumpFrac=0.1,offsetScale=3.0):
	kpno = coo.EarthLocation(lat=31.9583*u.degree,lon=-111.5967*u.degree,
	                         height=2120*u.m)
	utcOffset = -7*u.hour
	duration = duration*u.hour
	fixedElevation = fixedElevation*u.degree
	gbMin = 17*u.degree
	decMin = 30*u.degree
	decMax = 85*u.degree
	offsetScale *= u.degree
	#
	nAz = 1000
	altaz = coo.AltAz(alt=np.repeat(fixedElevation,nAz),
	                  az=np.linspace(0,360,nAz)*u.degree,
	                  obstime=startTime,location=kpno)
	gal = altaz.transform_to(coo.Galactic)
	cel = altaz.transform_to(coo.FK5)
	ii = np.where((gal.b > gbMin) & 
	              (cel.dec > decMin) & (cel.dec < decMax))[0]
	if len(ii)==0:
		return None
	np.random.shuffle(ii)
	i = ii[0]
	#
	az = altaz[i].az
	track = [cel[i]]
	times = [startTime]
	while times[-1]-startTime < duration:
		dAz = offsetScale*np.random.rand()
		if np.random.rand() < jumpFrac:
			dAz *= -3
		az = az + dAz
		t = times[-1] + shotTime * u.second
		a = coo.AltAz(alt=fixedElevation,az=az,obstime=t,location=kpno)
		c = a.transform_to(coo.FK5)
		g = a.transform_to(coo.Galactic)
		if c.dec < decMin or c.dec > decMax or g.b < gbMin:
			break
		track.append(c)
		times.append(t)
	return dict(coords=coo.SkyCoord(track),ut=times)

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
	minDuration = kwargs.get('minDuration',15*u.minute)
	utStart = kwargs.get('utStart','02:00')
	utEnd = kwargs.get('utEnd','13:00')
	utBreak = kwargs.get('utBreak',('05:00','10:00'))
	alltracks = {}
	shotTime = 100. #exposureTime + overheadTime
	for utdate,c in zip(utdates,itertools.cycle('rgbcymk')):
		alltracks[utdate] = []
		uttime = utStart
		utskip = True
		uts = utdate[:4]+'-'+utdate[4:6]+'-'+utdate[6:]
		uttime = Time(uts+' '+utStart)
		utend = Time(uts+' '+utEnd)
		utbreak = [Time(uts+' '+utBreak[0]),Time(uts+' '+utBreak[1])]
		nfail = 0
		while uttime < utend and nfail < 5:
			# random duration between 1/2 hour and 1.5 hours
			duration = 0.5 + np.random.rand()
			track = aztrack(uttime,shotTime,duration,elevation,**kwargs)
			if track is None:
				nfail += 1
				continue
			trackdur = track['ut'][-1]-track['ut'][0] 
			if trackdur < minDuration:
				# too short, try again
				nfail += 1
				continue
			nfail = 0
			plt.plot(track['coords'].ra.value,track['coords'].dec.value,
			         '%ss-'%c)
			uttime += trackdur
			alltracks[utdate].append(track)
			print uttime,trackdur.to(u.minute)
			if uttime > utbreak[0] and utskip:
				# skip a chunk in the middle of the night
				uttime = utbreak[1]
				utskip = False
	for ut in utdates:
		for j,t in enumerate(alltracks[ut]):
			print '%s %3d %7.2f %7.2f %7.2f' % \
			       (ut,j,t['coords'].dec[0].value,t['coords'].dec[-1].value,
			        (t['ut'][-1]-t['ut'][0]).to(u.minute))
	return alltracks

def formatut(ut,full=False):
	utd,utt = ut.iso.split()
	utd = utd.replace('-','')
	if full:
		return utd,utt[:-4]
	else:
		return utd,''.join(utt.split(':')[:2])

def dump_track(track,exposureTime,filt):
	from itertools import count
	utdate,utstr = formatut(track['ut'][0])
	scriptf = open('basscal_%s_%s_%s.txt' % (utdate,utstr,filt),'w')
	logf = open('basscal_%s_%s_%s.log' % (utdate,utstr,filt),'w')
	for i,ut,c in zip(count(),track['ut'],track['coords']):
		imtitle = 'cal%s%s%s_%02d' % (utdate,utstr,filt,i)
		# have to cut off last digit in dec for Bok
		cstr = c.to_string('hmsdms',sep='',precision=2)[:-1]
		scriptf.write("obs %.1f object '%s' 1 %s %s 2000.0\r\n" % 
		              (exposureTime,imtitle,filt,cstr))
		cstr = c.to_string('hmsdms',sep=':',precision=2)
		logf.write('%s   %s   %s\n' % (formatut(ut,True)[1],imtitle,cstr))
	scriptf.close()
	logf.close()

def make_track(*args,**kwargs):
	airmass = kwargs.get('airmass',1.4)
	exposureTime = kwargs.get('exposureTime',50.)
	overheadTime = kwargs.get('overheadTime',50.)
	duration = kwargs.get('duration',1.0)
	offsetScale = kwargs.get('offsetScale',3.0)
	jumpFrac = kwargs.get('jumpFrac',0.1)
	filt = kwargs.get('filt','g')
	if len(args)==1:
		t = Time.now()
		startTime = Time(t.iso.split()[0]+' '+args[0])
	elif len(args)==2:
		utdate,uttime = args
		utdate = utdate[:4]+'-'+utdate[4:6]+'-'+utdate[6:]
		startTime = Time(utdate+' '+uttime,format='iso')
	else:
		raise ValueError
	shotTime = exposureTime + overheadTime
	elevation = airmass2el(airmass)
	print
	print 'Constructing track at airmass=%.2f' % airmass
	print '  Duration: %.1f hr with exposure time = %.1fs (overhead %.1fs)' % \
	            (duration,exposureTime,overheadTime)
	print '  random params: offsetScale=%.1f deg, jumpFrac=%.1f' % \
	            (offsetScale,jumpFrac)
	print '  elevation is %.1f degrees' % elevation
	print '  start execution at ',startTime,' UT'
	print
	track = aztrack(startTime,shotTime,duration,
	                elevation,offsetScale=offsetScale,jumpFrac=jumpFrac)
	dump_track(track,exposureTime,filt)
	return track

if __name__=='__main__':
	import sys,getopt
	try:
		opts,args = getopt.getopt(sys.argv[1:],
		                          "a:d:e:o:",
		                         ['airmass=','duration=',
		                          'exposureTime=','filt=',
		                          'overheadTime=',])
	except getopt.GetoptError as err:
		print str(err)
		sys.exit(2)
	kwargs = {}
	for k,v in opts:
		kwargs[k.lstrip('--')] = float(v)
	t = make_track(*args,**kwargs)
#	ntry = 0
#	while ntry < 10:
#		t = make_track(*args,**kwargs)
#		if t is not None:
#			break
#		ntry += 1

