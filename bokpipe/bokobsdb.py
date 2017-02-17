#!/usr/bin/env python

import os,sys
import re
from glob import glob
from numpy import where
import fitsio
from math import cos,radians

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time,TimeDelta
from astropy.table import Table,vstack

badfloat = -9999.99

def generate_log(dirs,logFile,filters=None,objFilter=None,filePattern=None,
                 inTable=None,include_singlechip=False,
                 extraFields=(),extraTypes=(),extra_cb=None):
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
	t = Table(names=('frameIndex','utDir','fileName','utDate','date_obs',
	                 'imType','filter','objName','expTime',
	                 'ccdbin1','ccdbin2','oscan1','oscan2',
	                 'focusA','focusB','focusC',
	                 'hdrAirmass','airmass',
	                 'alt','az','ha','lst',
	                 'targetRa','targetDec','targetCoord',
	                 'utObs','mjdStart','mjdMid',
	                 'cameraTemp','dewarTemp',
	                 'outsideTemp','outsideHumidity','outsideDewpoint',
	                 'insideTemp','insideHumidity',
	                 'mirrorCellTemp','primaryTemp','strutTemp','primeTemp',
	                 'windSpeed','windDir','airTemp','relHumid','barom')
	                 +extraFields,
	          dtype=('i4','S15','S35','S8','S10',
	                 'S8','S8','S35','f4',
	                 'i4','i4','i4','i4',
	                 'f4','f4','f4',
	                 'f4','f4',
	                 'f4','f4','S10','S10',
	                 'f8','f4','S15',
	                 'S12','f8','f8',
	                 'f4','f4',
	                 'f4','f4','f4',
	                 'f4','f4',
	                 'f4','f4','f4','f4',
	                 'f4','f4','f4','f4','f4')
	                 +extraTypes)
	# enter the files as table rows
	for i,f in enumerate(files):
		fn = os.path.basename(f)
		fn = fn[:fn.find('.fits')]
		utDir = os.path.split(os.path.split(f)[0])[1]
		if inTable is not None:
			# assumes this combination is unique
			_i = where((inTable['fileName']==fn) &
			           (inTable['utDir']==utDir))[0]
			if len(_i)>0:
				print 'skipping ',f
				continue
		try:
			h = fitsio.read_header(f,0)
			h1 = fitsio.read_header(f,1)
		except:
			print 'ERROR: failed to open file: ',f
			continue
		if h['NCCDS']==1 and not include_singlechip:
			print 'WARNING: skipping single-chip image ',f
			continue
		row = []
		imageType = h.get('IMAGETYP','null').strip()
		filt = h.get('FILTER','null').strip()
		if filters is not None and filt not in filters:
			continue
		objName = str(h.get('OBJECT','null')).strip()
		if len(objName)==0:
			objName = 'null'
		if objFilter is not None and imageType == 'object' and \
		       not objFilter(objName):
			continue
		# airmass values in header are not very accurate (one decimal place)
		if not 'ELEVAT' in h:
			airmass = badfloat
		else:
			airmass = 1/cos(radians(90.-h['ELEVAT']))
			if 'AIRMASS' in h and abs(airmass-h['AIRMASS']) > 0.1:
				print airmass,h['AIRMASS'],h['ELEVAT']
				raise ValueError
		hdrAirmass = h.get('AIRMASS',badfloat)
		alt,az = h.get('ELEVAT',badfloat),h.get('AZIMUTH',badfloat)
		try:
			ha,lst = h['HA'].strip(),h['LST-OBS'].strip()
		except:
			ha,lst = '',''
		try:
			coord = h['RA'].rstrip()+' '+h['DEC'].rstrip()
			sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
			ra = sc.ra.degree
			dec = sc.dec.degree
		except:
			coord,ra,dec = '',badfloat,badfloat
		try:
			tObs = Time(h['DATE']+' '+h['UT'],scale='utc')
			# Add 5 hours and round down. Thus noon local (=7pm UT) becomes
			# midnight of the next day. This way afternoon cals get counted
			# with data from that night.
			tOff = tObs + TimeDelta(5*u.hour)
			# "2014-01-01 12:00:00.000" -> "20140101"
			utDate = str(tOff.utc).split()[0].replace('-','')
			mjdStart = tObs.mjd
			mjdMid = (tObs + TimeDelta(h['EXPTIME']*u.second)).mjd
		except:
			mjdStart,mjdMid = badfloat,badfloat
		try:
			focA,focB,focC = [float(x) for x in h['FOCUSVAL'].split('*')[1:]]
		except:
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
		if inTable is None:
			indx = i
		else:
			indx = len(inTable) + i
		row.extend([indx,utDir,fn,utDate,h.get('DATE-OBS','null')])
		row.extend([imageType,filt,objName,h.get('EXPTIME',-1)])
		row.extend([h.get('CCDBIN1',-1),h.get('CCDBIN2',-1),
		            h1.get('OVRSCAN1',-1),h1.get('OVRSCAN2',-1)])
		row.extend([focA,focB,focC])
		row.extend([hdrAirmass,airmass])
		row.extend([alt,az,ha,lst])
		row.extend([ra,dec,coord])
		row.extend([h.get('UT','null'),mjdStart,mjdMid])
		row.extend([h.get('CAMTEMP',-99999),h.get('DEWTEMP',-99999)])
		row.extend([outTemp,outHum,outDew])
		row.extend([inTemp,inHum])
		row.extend([mirrorCellTemp,primaryTemp,strutTemp,primeTemp])
		row.extend([windSpeed,windDir,airTemp,relHumid,barom])
		if len(extraFields)>0:
			row.extend([extra_cb(xf,h.get(xf)) for xf in extraFields])
		t.add_row(row)
		sys.stdout.write("\r%d/%d" % (i+1,len(files)))
		sys.stdout.flush()
	print
	if inTable is not None:
		t = vstack([inTable,t])
		t.sort('mjdStart')
	t.write(logFile,overwrite=True)

