#!/usr/bin/env python

#from __future__ import print_function

import os
import numpy as np
from astropy.table import Table,join,vstack

from bass import reform_filename,files2tiles,load_obsdb
import bokdepth

def get_bass_frames(obsLogFile):
	# read in the observations logs created from trawling the headers,
	# use the data common to the full frame as the starting point
	frames = Table.read(obsLogFile)
	# identify bass frames
	isbass = ( ((frames['filter']=='g')|(frames['filter']=='bokr')) &
	           (frames['imType']=='object')
	  )
	# XXX need to verify tile id!!!!
	#      this is a hacky way
	isbass &= [ len(s)==5 and s[-1] in '123' for s in frames['objName'] ]
	frames = frames[isbass]
	return frames

def process_dr3_frames(overwrite=False):
	import bokframechar
	frames = get_bass_frames('bassdr3frames.fits')
	# XXX use env
	rawDir = '/global/project/projectdirs/cosmo/staging/bok/BOK_Raw'
	dr3files = [ os.path.join(rawDir,f['utDir'],f['fileName']+'.fits.fz')
	                for f in frames ]
	outputdir = os.path.join(os.environ['SCRATCH'],'imageval','dr3')
	print 'processing ',len(dr3files),' images'
	bokframechar.process_images(dr3files,outputDir=outputdir,
	                            overwrite=overwrite)

#def match_results(ccds,etcVals):
#	ccds['seeing'][ii] = etcVals['seeing'][jj]
#	ccds['zpt'][ii] = etcVals['finalCal'][jj] + 25. \
#	                    - 2.5*np.log10(etcVals['expTime'][jj])
#	ccds['avsky'][ii] = etcVals['skyFlux'][jj]
#	ccds['ebv'][ii] = etcVals['E(B-V)'][jj]

def _temp_fn2expid_map(archivelistf):
	fnmap = {}
	with open(archivelistf) as archivef:
		for l in archivef:
			d = l.strip().split()
			noaofn = d[1]#.replace('.fz','').replace('.fits','')
			bokfn = d[2]
			try:
				expid = np.int32(bokfn[1:5]+bokfn[6:10])
			except:
				continue
			fnmap[noaofn] = expid
	return fnmap

def frames2ccds(frames,procdir,outfn='bass-ccds-annotated.fits',**kwargs):
	use_objname = kwargs.get('use_objname',False)
	framesOrig = frames
	frames = frames.copy()
	errlog = open(outfn.replace('.fits','_errors.log'),'w')
	# slice the columns that are needed
	cols = [ 'expTime', 'filter', 'date_obs', 'mjd', 'utObs',
	         'hdrAirmass', 'fileName', 'outsideHumidity', 'outsideTemp', ]
	frames = frames[cols]
	# translate some column names to match the DECaLS conventions
	_remap = [ ('expTime','exptime'), ('mjd','mjd_obs'), ('utObs','ut'),
	           ('hdrAirmass','airmass'), ('fileName','image_filename'), 
	           ('outsideHumidity','humidity'), ('outsideTemp','outtemp'),
	  ]
	for m in _remap:
		frames.rename_column(*m)
	# translate some data values
	frames['outtemp'] = 5.0/9.0 * (frames['outtemp']-32) # F->C
	frames['filter'][frames['filter']=='bokr'] = 'r'
	frames['propid'] = 'BASS'
	fns = np.core.defchararray.add(frames['image_filename'],'.fits.fz')
	if use_objname:
		# XXX this is really a hack for the Nov15 legacy fields...
		#     just undo their silly naming scheme!!!
		frames['image_filename'] = \
		   np.core.defchararray.add(framesOrig['objName'],'.fits')
		frames['image_filename'] = [_f.replace('bokr','r') for _f in frames['image_filename']]
	else:
		frames['image_filename'] = fns
	fnmap = _temp_fn2expid_map('nersc_noaoarchive_thru20160216.log') # XXX
	frames['expnum'] = [ fnmap[fn] for fn in fns ]
	frames['width'] = np.int32(4096)
	frames['height'] = np.int32(4032)
	# allocate dummy entries for per-ccd items
	frames['tileid'] = np.int32(0)
	frames['tilepass'] = np.uint8(0)
	for k in ['seeing','zpt','avsky','arawgain',
	          'crpix1','crpix2','cd1_1','cd1_2','cd2_1','cd2_2',
	          'ccdphrms','frameskyrms','ccdraoff','ccddecoff',
	          'ccdzpt','ccdmdncol',
	          'psfnorm_mean','ebv']:
		frames[k] = np.float32(0)
	for k in ['crval1','crval2']:
		frames[k] = np.float64(0)
	frames['ccdnum'] = np.int16(0)
	frames['ccdname'] = 'CCDX'
	frames['ccdnstar'] = np.int16(0)
	for amp in ['','a','b','c','d']:
		frames['ccdzpt'+amp] = np.float32(0)
		frames['ccdnmatch'+amp] = np.int16(0)
	# now explode the frames by 4 copies to make the per-CCD entries
	ccds = [ frames.copy() for i in range(4) ]
	for j,ccdNum in enumerate(range(1,5)):
		ccds[j]['ccdnum'] = ccdNum
		ccds[j]['ccdname'] = 'CCD%d' % ccdNum
	# values derived from images that are common to the whole frame
	perframekeys = ['seeing','zpt']
	# and those that need to be renamed from the meta-data file
	frmkmap = {'zpt':'zpim'}
	# ditto, but for values which are common to each CCD
	perccdkeys = ['arawgain','crpix1','crpix2','crval1','crval2',
	              'cd1_1','cd1_2','cd2_1','cd2_2',
	              'ccdzpt','ccdphrms','ccdnmatch','ccdmdncol',
	              'ccdraoff','ccddecoff']
	ccdkmap = {'ccdzpt':'zpccd','ccdphrms':'zprmsCCD',
	           'ccdnmatch':'nstarsccd','ccdmdncol':'medgicolor',
	           'ccdraoff':'raoff','ccddecoff':'decoff'}
	# ditto, but for values which are common to each amp
	perampkeys = ['ccdzpt','ccdnmatch']
	ampkmap = {'ccdzpt':'zpamp','ccdnmatch':'nstarsamp'}
	# iterate over the frames and extract values from the meta-data files
	for i,frame in enumerate(frames):
		fn = frame['image_filename']
		print 'processing ',fn
		# first extract the tile specifier from the log file
		objnm = framesOrig['objName'][i]
		try:
			tileid,tilepass = [int(objnm[:-1]),int(objnm[-1])]
			for ccd in ccds:
				ccd['tileid'][i] = tileid
				ccd['tilepass'][i] = tilepass
		except ValueError:
			# not a BASS tile
			pass
		#  now try to access processing results if they exist
		_fn = fn.replace('.fz','')
		metadatf = os.path.join(procdir,_fn.replace('.fits','.meta.npz'))
		try:
			metaDat = np.load(metadatf)
		except IOError:
			errlog.write('no processing for %s/%s\n' % 
			               tuple(framesOrig['utDir','fileName'][i]))
			continue
		for j,ccdNum in enumerate(range(1,5)):
			# even though these are per-frame still need a copy for each ccd
			ccds[j]['avsky'][i] = np.mean(metaDat['avsky'])
			for k in perframekeys:
				ccds[j][k][i] = metaDat[frmkmap.get(k,k)]
			for k in perccdkeys:
				ccds[j][k][i] = metaDat[ccdkmap.get(k,k)][j]
			for ai,aname in enumerate('abcd'):
				for k in perampkeys:
					ccds[j][k+aname][i] = metaDat[ampkmap.get(k,k)][j,ai]
		# finally extract PSF data
		psfexf = os.path.join(procdir,_fn.replace('.fits','.psf'))
		#   XXX a hacky workaround for file naming inconsistency
		if not os.path.exists(psfexf):
			psfexf += '.fits'
		try:
			psfs = bokdepth.make_PSFEx_psf_fromfile(psfexf,2048,2016)
		except IOError:
			errlog.write('no PSFEx model for %s/%s\n' % 
			               tuple(framesOrig['utDir','fileName'][i]))
			continue
		for j,ccdNum in enumerate(range(1,5)):
			psf = psfs[j]
			psf /= psf.sum()
			nea = np.sum(psf**2)**-1
			ccds[j]['psfnorm_mean'][i] = 1/nea
	# join all the CCDs into a single table, add a few more fields
	allccds = vstack(ccds)
	# the WCS used sets the center of the focal plane in crvaln, so this
	# works, but watch out if it changes...
	allccds['ra_bore'] = allccds['crval1']
	allccds['dec_bore'] = allccds['crval2']
	# finally sort by exposure number and then ccd number
	allccds.sort(['expnum','ccdnum'])
	allccds.write(outfn,overwrite=True)
	errlog.close()

def _tmp_dr3():
	#frames = Table.read('../basschute/config/nov2015_mod.fits')
	#  XXX need to update this by copying in "good" field from _mod
	frames = Table.read('../basschute/config/nov2015_new.fits')
	ii = np.concatenate([np.arange(33,70),np.arange(98,111)])
	frames = frames[ii]
	frames = frames[frames['objName'] != 'null']
	outf = 'bass-ccds-idmnov2015.fits'
	frames2ccds(frames,'/global/homes/i/imcgreer/bok/reduced/nov15data_ian/',
	            outf,use_objname=True)
	t = Table.read(outf)
	for s in ['','a','b','c','d']:
		t['ccdzpt'+s] += -2.5*np.log10(t['arawgain'])
	t.write(outf,overwrite=True)

if __name__=='__main__':
	import sys
	if len(sys.argv)>1:
		procdir = os.environ['BASSFRAMEDIR']
		frames = get_bass_frames(sys.argv[1])
		frames2ccds(frames,procdir)
	else:
		process_dr3_frames()

