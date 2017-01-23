#!/usr/bin/env python

#from __future__ import print_function

import os
import time
import numpy as np
from astropy.table import Table,join,vstack
from astropy.io import fits

from bass import reform_filename,files2tiles,load_obsdb
import bokdepth
from bokframechar import get_ps1_zpastrom

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
	#isbass &= [ len(s)==5 and s[-1] in '123' for s in frames['objName'] ]
	isbass &= (frames['expTime'] > 20)
	frames = frames[isbass]
	return frames

def process_frames(frames,overwrite=False,nproc=5,searchrad=None,dr3=False,
                   rawDir=None):
	import multiprocessing
	import bokframechar
	if rawDir is None:
		rawDir = os.environ['BASSRAWDATA']
	files = [ os.path.join(rawDir,f['utDir'],f['fileName']+'.fits.fz')
	                for f in frames ]
	outputdir = os.path.join(os.environ['SCRATCH'],'imageval')
	if dr3:
		outputdir = os.path.join(outputdir,'dr3')
	else:
		outputdir = os.path.join(outputdir,'test')
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	print 'processing ',len(files),' images'
	kwargs = {'outputDir':outputdir,'overwrite':overwrite}
	if searchrad:
		kwargs['searchrad'] = searchrad
	if nproc==1:
		bokframechar.process_raw_images(files,**kwargs)
	else:
		for i in range(nproc):
			_files = files[i::nproc]
			p = multiprocessing.Process(target=bokframechar.process_raw_images,
			                            args=(_files,),kwargs=kwargs)
			p.start()

def _extract_metadata_fromheader(ccds,i,expnum,oframe,rdxDir=None):
	expstr = str(expnum)
	if rdxDir is None:
		rdxDir = os.path.join(os.environ['BASSDATA'],'reduced')
	fdir = oframe['utDir']
	# XXX shouldn't be here...
	# gains taken from Mike Lesser's doc on 23-Sep-15
	gains = [ [1.39,1.72,1.31,1.46],
	          [1.48,1.43,1.53,1.48],
	          [1.37,1.35,1.43,1.53],
	          [1.35,1.31,1.66,1.39] ]
###	if not os.path.exists(os.path.join(rdxdir,fdir)):
###		# some of the dirs have been appended with the filter
###		fdir += oframe['filter'][-1]
	for ccdNum,ccd in enumerate(ccds,start=1):
		filt = oframe['filter'].replace('bokr','r')
		imfn = ''.join(['p',expstr[:4],filt,expstr[4:],
		                '_%d'%ccdNum,'.fits'])
		fpath = os.path.join(fdir,imfn)
		# extensions are not always named so use number
		hdr = fits.getheader(os.path.join(rdxDir,fpath),0)#'CCD%d'%ccdNum)
		for k in ['seeing','crpix1','crpix2','crval1','crval2',
		          'cd1_1','cd1_2','cd2_1','cd2_2','cali_ref']:
			ccd[k][i] = hdr[k.upper()]
		for k in ['zpt','zpta','zptb','zptc','zptd','phrms','nstar']:
			ccd['ccd'+k][i] = hdr['CCD'+k.upper()]
		for k in ['nmatch','nmatcha','nmatchb','nmatchc','nmatchd',
		          'mdncol','skyrms','raoff','decoff']:
			ccd['ccd'+k][i] = hdr[k.upper()]
		try:
			ccd['avsky'][i] = hdr['SKADU']
		except KeyError:
			pass
		ccd['arawgain'][i] = np.mean(gains[ccdNum-1])
		ccd['image_filename'][i] = fpath
	return ccds

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

def _extract_metadata_fromfile(ccds,i,metadatf):
	metaDat = np.load(metadatf)
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
	return ccds

def _get_naoc_name(expnum,ccdNum,oframe,procdir,sfx=''):
	expstr = str(expnum)
	filt = oframe['filter'].replace('bokr','r')
	imfn = ''.join(['p',expstr[:4],filt,expstr[4:],'_%d'%ccdNum,
	                sfx,'.fits'])
	return os.path.join(procdir,oframe['utDir'],imfn)

def _redo_zeropoints(ccds,i,expnum,oframe,procdir,rdxDir=None):
	if rdxDir is None:
		rdxDir = os.path.join(os.environ['BASSDATA'],'reduced')
	for j,ccdNum in enumerate(range(1,5)):
		filt = oframe['filter'].replace('bokr','r')
		#expTime = oframe['expTime']
		ps1mf = _get_naoc_name(expnum,ccdNum,oframe,procdir,'.ps1match')
		ps1m = fits.getdata(ps1mf)
		# set expTime=1.0 because image units are e-/s
		zps = get_ps1_zpastrom(ps1m,filt,1.0)#expTime)
		# even though these are per-frame still need a copy for each ccd
		for k in perframekeys:
			if k=='seeing':
				# since this is currently calculated per-ccd and only one ccd
				# is available here, not actually using the mean across ccds
				ccds[j][k][i] = zps[frmkmap.get(k,k)][j]
			else:
				ccds[j][k][i] = zps[frmkmap.get(k,k)]
		# same as for seeing
		ccds[j]['avsky'][i] = zps['avsky'][j]
		for k in perccdkeys:
			if k.startswith('ccd'):
				ccds[j][k][i] = zps[ccdkmap.get(k,k)][j]
		for ai,aname in enumerate('abcd'):
			for k in perampkeys:
				ccds[j][k+aname][i] = zps[ampkmap.get(k,k)][j,ai]
	return ccds

def frames2ccds(frames,procdir,outfn='bass-ccds-annotated.fits',**kwargs):
	imgsource = kwargs.get('imgsource','raw')
	redozeropts = kwargs.get('redozeropts',False)
	rdxDir = kwargs.get('reduxdir',
	                    os.path.join(os.environ['BASSDATA'],'reduced'))
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
	frames['camera'] = '90prime'
	fns = np.core.defchararray.add(frames['image_filename'],'.fits.fz')
	if imgsource=='nov2015idm':
		# XXX this is really a hack for the Nov15 legacy fields...
		#     just undo their silly naming scheme!!!
		frames['image_filename'] = \
		   np.core.defchararray.add(framesOrig['objName'],'.fits')
		frames['image_filename'] = [_f.replace('bokr','r') for _f in frames['image_filename']]
	elif imgsource=='reduced':
		pass
	elif imgsource=='raw':
		frames['image_filename'] = fns
	else:
		raise ValueError('imgsource %s not recognized' % imgsource)
	frames['expnum'] = [ np.int32(fn[1:5]+fn[6:10]) 
	                           for fn in framesOrig['DTACQNAM'] ]
	frames['width'] = np.int32(4096)
	frames['height'] = np.int32(4032)
	# allocate dummy entries for per-ccd items
	frames['tileid'] = np.int32(0)
	frames['tilepass'] = np.uint8(0)
	for k in ['seeing','zpt','avsky','arawgain',
	          'crpix1','crpix2','cd1_1','cd1_2','cd2_1','cd2_2',
	          'ccdphrms','ccdskyrms','ccdraoff','ccddecoff',
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
	frames['cali_ref'] = 'NULLCAL'
	# now explode the frames by 4 copies to make the per-CCD entries
	ccds = [ frames.copy() for i in range(4) ]
	for j,ccdNum in enumerate(range(1,5)):
		ccds[j]['ccdnum'] = ccdNum
		ccds[j]['ccdname'] = 'CCD%d' % ccdNum
	# iterate over the frames and extract values from the meta-data files
	for i,frame in enumerate(frames):
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
		print 'processing frame ',i+1,frame['expnum'],frame['image_filename']
		# now try to extract processing results if they exist
		try:
			if imgsource=='reduced':
				ccds = _extract_metadata_fromheader(ccds,i,frame['expnum'],
				                                    framesOrig[i],
				                                    rdxDir=rdxDir)
			else:
				fn = frame['image_filename']
				_fn = fn.replace('.fz','')
				metadatf = os.path.join(procdir,
				                        _fn.replace('.fits','.meta.npz'))
				ccds = _extract_metadata_fromfile(ccds,i,metadatf)
		except (IOError,KeyError) as e:
			errlog.write('no processing for %s/%s\n' % 
			               tuple(framesOrig['utDir','fileName'][i]))
			continue
		if redozeropts:
			try:
				ccds = _redo_zeropoints(ccds,i,frame['expnum'],
				                        framesOrig[i],procdir,rdxDir=rdxDir)
			except (IOError) as e:
				errlog.write('no PS1 match data for %s/%s\n' % 
				               tuple(framesOrig['utDir','fileName'][i]))
				continue
		# finally extract PSF data
		try:
			if imgsource=='reduced':
				if redozeropts:
					# very hacky way to access regenerated PSF files
					psffns = [ _get_naoc_name(frame['expnum'],ccdNum,
					                          framesOrig[i],procdir,
					                          '.psf').replace('.fits','')
					            for ccdNum,ccd in enumerate(ccds,start=1) ]
				else:
					psffns = [ os.path.join(rdxDir,
					       ccd['image_filename'][i].replace('.fits','.psf'))
					            for ccd in ccds ]
				psfs = [ bokdepth.make_PSFEx_psf_fromfile(psffn,2048,2016)
				            for psffn in psffns ]
			else:
				psfexf = os.path.join(procdir,_fn.replace('.fits','.psf'))
				#   XXX a hacky workaround for file naming inconsistency
				if not os.path.exists(psfexf):
					psfexf += '.fits'
				psfs = bokdepth.make_PSFEx_psf_fromfile(psfexf,2048,2016)
		except:
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
	if imgsource=='reduced':
		allccds['image_hdu'] = 0
	else:
		allccds['image_hdu'] = allccds['ccdnum']
	# finally sort by exposure number and then ccd number
	allccds.sort(['expnum','ccdnum'])
	allccds.write(outfn,overwrite=True)
	errlog.close()

def _tmp_dr3():
	#frames = Table.read('../basschute/config/nov2015_mod.fits')
	#  XXX need to update this by copying in "good" field from _mod
	#frames = Table.read('../basschute/config/nov2015_new.fits')
	frames = Table.read('bass_Nov2015toFeb2016.fits')
	ii = np.concatenate([np.arange(33,70),np.arange(98,111)])
	frames = frames[ii]
	frames = frames[frames['objName'] != 'null']
	outf = 'bass-ccds-idmnov2015.fits'
	frames2ccds(frames,'/global/homes/i/imcgreer/bok/reduced/nov15data_ian/',
	            outf,imgsource='nov2015idm')
	t = Table.read(outf)
	for s in ['','a','b','c','d']:
		t['ccdzpt'+s] += -2.5*np.log10(t['arawgain'])
	t.write(outf,overwrite=True)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("inputFiles",type=str,nargs='*',
	       help="input FITS tables with image entries")
	parser.add_argument("-f","--ccdsfile",type=str,
	                    default='bass-ccds-tmp.fits',
	       help="name of output CCDs file")
	parser.add_argument("-n","--framenums",type=str,
	       help="subset of frames to process (e.g., '1,200')")
	parser.add_argument("-o","--outputdir",type=str,
	                    default=os.environ.get('BASSFRAMEDIR'),
	       help="directory containing processed data [default:$BASSFRAMEDIR]")
	parser.add_argument("-p","--process",action="store_true",
	       help="run processing on raw images instead of making CCDs file")
	parser.add_argument("--check",action="store_true",
	       help="check processing status instead of making CCDs file")
	parser.add_argument("-r","--raw",action="store_true",
	                    help="input is raw images, not processed")
	parser.add_argument("-R","--redo",action="store_true",
	       help="reprocess and overwrite existing files")
	parser.add_argument("-m","--multiproc",type=int,default=1,
	       help="number of processes to launch if doing processing")
	parser.add_argument("-v","--verbose",action="store_true",
	       help="increase verbosity")
	parser.add_argument("-x","--exclude",type=str,
	       help="file with list of images to exclude")
	parser.add_argument("--dr3",action="store_true",
	       help="select only DR3 frames")
	parser.add_argument("--searchradius",type=float,
	       help="WCS search radius in deg when processing")
	parser.add_argument("--redozeropts",action="store_true",
	                    help="redo PS1 zeropoints")
	parser.add_argument("--reduxdir",type=str,
	       help="location of processed data")
	args = parser.parse_args()
	frames = vstack([get_bass_frames(t) for t in args.inputFiles])
	if args.framenums is not None:
		i1,i2 = [int(v) for v in args.framenums.split(',')]
		frames = frames[i1-1:i2]
	if args.exclude is not None:
		with open(args.exclude) as xf:
			badlist = { l.strip():True for l in xf }
		drop = []
		for i,f in enumerate(frames):
			if f['fileName'] in badlist:
				drop.append(i)
		frames.remove_rows(drop)
	if args.process:
		kwargs = {'searchrad':args.searchradius}
		process_frames(frames,overwrite=args.redo,
		               nproc=args.multiproc,dr3=args.dr3,**kwargs)
	else:
		imgsrc = 'raw' if args.raw else 'reduced'
		if args.check:
			nfound = 1
			for i,_f in enumerate(frames,start=1):
				_fn = _f['fileName']
				metadatf = os.path.join(args.outputdir,_fn+'.meta.npz')
				if os.path.exists(metadatf):
					if args.verbose:
						t = time.ctime(os.path.getmtime(metadatf))
						print '%4d %s %s' % (i,_fn,t)
					nfound += 1
			print nfound,len(frames)
		else:
			frames2ccds(frames,args.outputdir,args.ccdsfile,
			            imgsource=imgsrc,redozeropts=args.redozeropts,
			            reduxdir=args.reduxdir)

