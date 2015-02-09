#!/usr/bin/env python

import os
import itertools
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.time import Time
import astropy.coordinates as coo
from astropy import units as u

import matplotlib.pyplot as plt

utcOffset = -7*u.hour
decMin = 30*u.degree
decMax = 85*u.degree

def get_desi_tiles():
	tiles = fits.getdata('desi-tiles.fits')
	ii = np.where((tiles.IN_DESI > 0)  & (tiles.DEC > 30))[0]
	return tiles[ii]

def get_desi_bounds():
	tiles = get_desi_tiles()
	stripewidth = 10.0
	pad = 5.0
	decedges = np.arange(25,95,stripewidth)
	ii = np.digitize(tiles.DEC,decedges)
	nbins = len(decedges)-1
	ramin = [tiles.RA[ii==i+1].min() for i in range(nbins)]
	ramax = [tiles.RA[ii==i+1].max() for i in range(nbins)]
	ramin = np.array(ramin+[120]) - pad
	ramax = np.array(ramax+[250]) + pad
	decbins = decedges + stripewidth/2
	#
	dec1 = np.arange(30,86.5,1.0) * u.degree
	ra1 = np.empty_like(dec1)
	ra2 = np.empty_like(dec1)
	for i in range(len(dec1)):
		j = int((dec1[i] - decMin).value//stripewidth)
		y0,y1 = decbins[[j,j+1]] * u.degree
		x0,x1 = ramin[[j,j+1]] * u.degree
		ra1[i] = x0 + ((x1-x0)/(y1-y0))*(dec1[i]-y0)
		x0,x1 = ramax[[j,j+1]] * u.degree
		ra2[i] = x0 + ((x1-x0)/(y1-y0))*(dec1[i]-y0)
	return dict(decStripe=dec1,raMin=ra1,raMax=ra2)

def in_desi_bounds(ra,dec,bounds):
	try:
		ii = np.digitize(dec,bounds['decStripe'])
	except ValueError:
		ii = np.where((dec > bounds['decStripe'][:-1]) & 
		              (dec < bounds['decStripe'][1:]))[0] + 1
	return ( (ra > bounds['raMin'][ii-1]) & (ra < bounds['raMax'][ii-1]) &
	         (dec > decMin) & (dec < decMax) )

def get_kpno():
	return coo.EarthLocation(lat=31.9583*u.degree,lon=-111.5967*u.degree,
	                         height=2120*u.m)

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
	kpno = get_kpno()
	duration = duration*u.hour
	fixedElevation = fixedElevation*u.degree
	offsetScale *= u.degree
	# find positions within the azimuth track at the given elevation
	# that are within the survey footprint
	nAz = 100
	altaz = coo.AltAz(alt=np.repeat(fixedElevation,nAz),
	                  az=np.linspace(0,360,nAz)*u.degree,
	                  obstime=startTime,location=kpno)
	cel = altaz.transform_to(coo.FK5)
	bounds = get_desi_bounds()
	ii = np.where(in_desi_bounds(cel.ra,cel.dec,bounds))[0]
	if len(ii)==0:
		return None
	# pick the starting point randomly from the list of valid positions
	np.random.shuffle(ii)
	i = ii[0]
	# now build up the track by making random azimuthal offsets
	az = altaz[i].az
	coords = [cel[i]]
	azimuths = [az]
	times = [startTime]
	while times[-1]-startTime < duration:
		dAz = offsetScale*np.random.rand()
		if np.random.rand() < jumpFrac:
			dAz *= -3
		az = az + dAz
		t = times[-1] + shotTime * u.second
		a = coo.AltAz(alt=fixedElevation,az=az,obstime=t,location=kpno)
		c = a.transform_to(coo.FK5)
		# stop if this move takes the track out of bounds
		if not in_desi_bounds(c.ra,c.dec,bounds):
			break
		coords.append(c)
		azimuths.append(az)
		times.append(t)
	return dict(coords=coo.SkyCoord(coords),ut=times,azimuths=azimuths)

def airmass2el(airmass):
	return 90 - np.degrees(np.arccos(1/airmass))

def plot_gal_bound(gb=17,c='r'):
	gl = np.linspace(-180,180,100) * u.degree
	gb = np.repeat(gb,len(gl)) *u.degree
	gbound = coo.SkyCoord(gl,gb,coo.Galactic)
	bound = gbound.transform_to(coo.FK5)
	plt.plot(bound.ra,bound.dec,c=c)

def plot_desi_tiles():
	tiles = get_desi_tiles()
	plt.scatter(tiles.RA,tiles.DEC,c='0.2',marker='+',s=10)

def desi_airmasses():
	tiles = get_desi_tiles()
	utStart = '02:00'
	utEnd = '13:00'
	dt = 0.2 * u.hour
	kpno = get_kpno()
	ntiles = len(tiles)
	tiles = coo.SkyCoord(tiles.RA*u.degree,tiles.DEC*u.degree,coo.FK5)
	for month in range(2,8):
		tStart = Time('2015-%02d-03 %s' % (month,utStart))
		tEnd = Time('2015-%02d-03 %s' % (month,utEnd))
		t = tStart
		while t<tEnd:
			bokframe = coo.AltAz(obstime=t,location=kpno)
			a = tiles.transform_to(bokframe)
			print a.az.value
			print t
			t += dt
		break

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
	if os.path.exists('desi-tiles.fits'):
		plot_desi_tiles()
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
			         c=c,marker='o',mfc='none',ms=4,mec=c)
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
			        (t['ut'][-1]-t['ut'][0]).to(u.minute).value)
	return alltracks

def formatut(ut,full=False):
	utd,utt = ut.iso.split()
	utd = utd.replace('-','')
	if full:
		return utd,utt[:-4]
	else:
		return utd,''.join(utt.split(':')[:2])

def dump_track(track,exposureTime,filt,nowrite):
	from itertools import count
	utdate,utstr = formatut(track['ut'][0])
	if not nowrite:
		scriptf = open('basscal_%s_%s_%s.txt' % (utdate,utstr,filt),'w')
		logf = open('basscal_%s_%s_%s.log' % (utdate,utstr,filt),'w')
	for i,ut,c,az in zip(count(),track['ut'],track['coords'],track['azimuths']):
		imtitle = 'cal%s%s%s_%02d' % (utdate,utstr,filt,i)
		# have to cut off last digit in dec for Bok
		cstr = c.to_string('hmsdms',sep='',precision=2)[:-1]
		if not nowrite:
			scriptf.write("obs %.1f object '%s' 1 %s %s 2000.0\r\n" % 
			              (exposureTime,imtitle,filt,cstr))
		cstr = c.to_string('hmsdms',sep=':',precision=2)
		if not nowrite:
			logf.write('%s MST  %s UT   %s   %s   %5.2f\n' % 
			           (formatut(ut+utcOffset,True)[1],formatut(ut,True)[1],
			            imtitle,cstr,az.value))
		print '%s MST   %s UT   %s   %s   %5.2f' %  \
		           (formatut(ut+utcOffset,True)[1],formatut(ut,True)[1],
		            imtitle,cstr,az.value)
	if not nowrite:
		scriptf.close()
		logf.close()

def make_track(startTime,**kwargs):
	airmass = kwargs.get('airmass',1.4)
	exposureTime = kwargs.get('exposureTime',50.)
	overheadTime = kwargs.get('overheadTime',50.)
	duration = kwargs.get('duration',1.0)
	offsetScale = kwargs.get('offsetScale',3.0)
	jumpFrac = kwargs.get('jumpFrac',0.1)
	filt = kwargs.get('filt','g')
	shotTime = exposureTime + overheadTime
	elevation = airmass2el(airmass)
	print
	print 'Constructing track at airmass=%.2f' % airmass
	print '  Duration: %.1f hr with exposure time = %.1fs (overhead %.1fs)' % \
	            (duration,exposureTime,overheadTime)
	print '  random params: offsetScale=%.1f deg, jumpFrac=%.1f' % \
	            (offsetScale,jumpFrac)
	print '  elevation is %.1f degrees' % elevation
	print '  start execution at ',startTime+utcOffset,' MST'
	print
	track = aztrack(startTime,shotTime,duration,
	                elevation,offsetScale=offsetScale,jumpFrac=jumpFrac)
	if track is not None:
		dump_track(track,exposureTime,filt,kwargs.get('nowrite',False))
	return track

if __name__=='__main__':
	import sys,argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('utdate',nargs='?',default=None,
	                    help="UT date of observation, in format YYYYMMDD "+
	                         "[default: today's date]")
	parser.add_argument('uttime',#nargs='?',default='+5',
	              help="UT time of observation, in format HH:MM[:SS] "+
	                   "or offset from current time in minutes, e.g., "+
	                   "'+5' mean 5 minutes from now")# [this is the default]")
	parser.add_argument("-a","--airmass",type=float,default=1.4,
	              help="airmass of observation [default=1.4]")
	parser.add_argument("-d","--duration",type=float,default=1.0,
	              help="duration of observation in hours [default=1.0]")
	parser.add_argument("-e","--exposureTime",type=float,default=50,
	              help="exposure time in seconds [default=50]")
	parser.add_argument("-f","--filt",default='g',
	              help="filter name [default=g]")
	parser.add_argument("-o","--overheadTime",type=float,default=50,
	              help="overhead time in seconds [default=50]")
	parser.add_argument("--offsetScale",type=float,default=3.0,
	              help="scale of random offsets in degrees [default=3]")
	parser.add_argument("--jumpFrac",type=float,default=0.1,
	              help="fraction of offsets to be negative in azimuth "+
	                   "[default=0.1]")
	parser.add_argument("--nowrite",action="store_true",
	              help="skip writing output files, just dump to screen")
	args = parser.parse_args()
	# determine the starting time for the sequence
	t = Time.now()
	if args.utdate is None:
		utdate = t.iso.split()[0]
	else:
		utdate = args.utdate[:4]+'-'+args.utdate[4:6]+'-'+args.utdate[6:]
	if args.uttime.startswith('+'):
		dt = float(args.uttime[1:]) * u.minute
		print 'dt = ',dt
		startTime = t + dt
	else:
		startTime = Time(utdate+' '+args.uttime,format='iso')
	# convert the optional arguments to a dictionary and pass to make_track
	opts = vars(args)
	ignore = opts.pop('utdate')
	ignore = opts.pop('uttime')
	kwargs = { k : opts[k] for k in opts if opts[k] != None }
	# loop a few times in case the random starting point results in a failed
	# sequence (i.e., goes out of bounds), or if the track fails otherwise
	# (e.g., if no part of the footprint is visible at the specified airmass
	#  for the specified observing time)
	# [it would be better to recognize the latter case and give up]
	ntry = 0
	while ntry < 10:
		t = make_track(startTime,**kwargs)
		if t is not None:
			break
		ntry += 1

