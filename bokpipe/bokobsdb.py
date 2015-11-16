#!/usr/bin/env python

import os,sys
import re
from glob import glob
import fitsio
from math import cos,radians

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table

badfloat = -9999.99

def generate_log(dirs,logFile,filters=None,objFilter=None,filePattern=None):
	# load all FITS files in the specified directories
	if filePattern is None:
		filePattern = '*.fits*'
	if type(dirs) is str:
		files = glob(os.path.join(dirs,filePattern))
	elif type(dirs) is list:
		files = []
		for d in dirs:
			files.extend(glob(os.path.join(d,filePattern)))
	files.sort()
	# create the storage table
	t = Table(names=('frameIndex','utDir','fileName','utDate',
	                 'imType','filter','objName','expTime',
	                 'ccdbin1','ccdbin2','focusA','focusB','focusC',
	                 'hdrAirmass','airmass',
	                 'alt','az','ha','lst',
	                 'targetRa','targetDec','targetCoord',
	                 'utObs','mjd',
	                 'cameraTemp','dewarTemp',
	                 'outsideTemp','outsideHumidity','outsideDewpoint',
	                 'insideTemp','insideHumidity',
	                 'mirrorCellTemp','primaryTemp','strutTemp','primeTemp',
	                 'windSpeed','windDir','airTemp','relHumid','barom'),
	          dtype=('i4','S15','S35','S8',
	                 'S8','S8','S15','f4',
	                 'i4','i4','f4','f4','f4',
	                 'f4','f4',
	                 'f4','f4','S10','S10',
	                 'f8','f4','S15',
	                 'S12','f8',
	                 'f4','f4',
	                 'f4','f4','f4',
	                 'f4','f4',
	                 'f4','f4','f4','f4',
	                 'f4','f4','f4','f4','f'))
	# enter the files as table rows
	for i,f in enumerate(files):
		try:
			h = fitsio.read_header(f,0)
		except:
			print 'ERROR: failed to open file: ',f
			continue
		row = []
		fn = os.path.basename(f).rstrip('.gz').rstrip('.fits')
		utDir = os.path.split(os.path.split(f)[0])[1]
		imageType = h['IMAGETYP'].strip()
		filt = h['FILTER'].strip()
		if filters is not None and filt not in filters:
			continue
		objName = h['OBJECT'].strip()
		if len(objName)==0:
			objName = 'null'
		if objFilter is not None and imageType == 'object' and \
		       not objFilter(objName):
			continue
		utDate = h['DATE'].strip().replace('-','')
		# airmass values in header are not very accurate (one decimal place)
		if not 'ELEVAT' in h:
			airmass = badfloat
		else:
			airmass = 1/cos(radians(90.-h['ELEVAT']))
			if abs(airmass-h['AIRMASS']) > 0.1:
				print airmass,h['AIRMASS'],h['ELEVAT']
				raise ValueError
		try:
			hdrAirmass = h['AIRMASS']
		except ValueError:
			hdrAirmass = badfloat
		try:
			alt,az = h['ELEVAT'],h['AZIMUTH'] 
		except ValueError:
			alt,az = badfloat,badfloat
		try:
			ha,lst = h['HA'].strip(),h['LST-OBS'].strip()
		except ValueError:
			ha,lst = '',''
		try:
			coord = h['RA'].rstrip()+' '+h['DEC'].rstrip()
			sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
			ra = sc.ra.degree
			dec = sc.dec.degree
		except ValueError:
			coord,ra,dec = '',badfloat,badfloat
		try:
			mjd = Time(h['DATE']+' '+h['UT'],scale='utc').mjd
		except ValueError:
			mjd = badfloat
		try:
			focA,focB,focC = [float(x) for x in h['FOCUSVAL'].split('*')[1:]]
		except ValueError:
			focA,focB,focC = badfloat,badfloat,badfloat
		try:
			tempstr = r'.*TEMP_F=(.*) HUMID_%=(.*) DEWPOINT_F=(.*)'
			m = re.match(tempstr,h['TEMPS0'])
			outTemp,outHum,outDew = [float(x) for x in m.groups()]
		except:
			outTemp,outHum,outDew = badfloat,badfloat,badfloat
		try:
			m = re.match(tempstr,h['TEMPS1'])
			inTemp,inHum,inDew = [float(x) for x in m.groups()]
		except:
			inTemp,inHum,inDew = badfloat,badfloat,badfloat
		try:
			m = re.match(tempstr,h['TEMPS2'])
			mirrorCellTemp,mcHum,mcDew = [float(x) for x in m.groups()]
		except:
			mirrorCellTemp,mcHum,mcDew = badfloat,badfloat,badfloat
		try:
			tempstr = r'.*TEMP_F=(.*)'
			m = re.match(tempstr,h['TEMPS4'])
			primaryTemp = float(m.groups()[0])
		except:
			primaryTemp = badfloat
		try:
			m = re.match(tempstr,h['TEMPS5'])
			strutTemp = float(m.groups()[0])
		except:
			strutTemp = badfloat
		try:
			m = re.match(tempstr,h['TEMPS6'])
			primeTemp = float(m.groups()[0])
		except:
			primeTemp = badfloat
		try:
			wstr = r'wind_speed=(.*)'
			m = re.match(wstr,h['WEATHER0'])
			windSpeed = float(m.groups()[0])
		except:
			windSpeed = badfloat
		try:
			wstr = r'wind_direction=(.*)'
			m = re.match(wstr,h['WEATHER1'])
			windDir = float(m.groups()[0])
		except:
			windDir = badfloat
		try:
			wstr = r'air_temperature=(.*)'
			m = re.match(wstr,h['WEATHER2'])
			airTemp = float(m.groups()[0])
		except:
			airTemp = badfloat
		try:
			wstr = r'relative_humid=(.*)'
			m = re.match(wstr,h['WEATHER3'])
			relHumid = float(m.groups()[0])
		except:
			relHumid = badfloat
		try:
			wstr = r'barometer=(.*) in'
			m = re.match(wstr,h['WEATHER4'])
			barom = float(m.groups()[0])
		except:
			barom = badfloat
		#
		row.extend([i,utDir,fn,utDate])
		row.extend([imageType,filt,objName,h['EXPTIME']])
		row.extend([h['CCDBIN1'],h['CCDBIN2'],focA,focB,focC])
		row.extend([hdrAirmass,airmass])
		row.extend([alt,az,ha,lst])
		row.extend([ra,dec,coord])
		row.extend([h['UT'],mjd])
		row.extend([h['CAMTEMP'],h['DEWTEMP']])
		row.extend([outTemp,outHum,outDew])
		row.extend([inTemp,inHum])
		row.extend([mirrorCellTemp,primaryTemp,strutTemp,primeTemp])
		row.extend([windSpeed,windDir,airTemp,relHumid,barom])
		t.add_row(row)
		sys.stdout.write("\r%d/%d" % (i+1,len(files)))
		sys.stdout.flush()
	print
	if os.path.exists(logFile):
		os.unlink(logFile)
	t.write(logFile)

