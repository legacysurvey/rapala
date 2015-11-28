#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.stats import sigma_clip

import bokrmpipe
from bokrmgnostic import srcor

def aperture_phot(file_map,inputType='sky',**kwargs):
	from astropy.table import Table,vstack
	redo = kwargs.get('redo',False)
	aperRad = np.concatenate([np.arange(2,9.51,1.5),[15.,22.5]])
	sdss = fitsio.read(os.environ['BOK90PRIMEDIR']+'/../data/sdss.fits',1)
	bpMask = file_map('MasterBadPixMask4')
	catDir = os.path.join(file_map.procDir,'catalogs')
	if not os.path.exists(catDir):
		os.mkdir(catDir)
	catPfx = 'bokrm_sdss'
	for filt in file_map.iterFilters():
		for utd in file_map.iterUtDates():
			fn = '.'.join([catPfx,utd,filt,'cat','fits'])
			catFile = os.path.join(catDir,fn)
			if os.path.exists(catFile):
				if redo:
					os.unlink(catFile)
				else:
					print catFile,' already exists, skipping'
					continue
			files,frames = file_map.getFiles(imType='object',
			                                 with_frames=True)
			if files is None:
				continue
			allPhot = []
			for imFile,frame in zip(files,frames):
				imageFile = file_map(inputType)(imFile)
				aHeadFile = imageFile.replace('.fits','.ahead')
				print 'processing ',imageFile
				phot = bokphot.aper_phot_image(imageFile,
				                               sdss['ra'],sdss['dec'],
				                               aperRad,bpMask,
				                               aHeadFile=aHeadFile,**kwargs)
				phot['frameNum'] = np.int32(frame)
				if phot is None:
					print 'no apertures found!!!!'
					continue
				allPhot.append(phot)
			allPhot = vstack(allPhot)
			allPhot.write(catFile)

def calc_zero_points(dataMap,magRange=(16.,19.5),aperNum=-2):
	pfx = 'bokrm_sdss'
	aperCatDir = os.path.join(dataMap.procDir,'catalogs')
	sdss = fits.getdata(os.environ['BOK90PRIMEDIR']+'/../data/sdss.fits',1)
	for filt in dataMap.iterFilters():
		is_mag = ( (sdss[filt]>=magRange[0]) & (sdss[filt]<=magRange[1]) )
		ref_ii = np.where(is_mag)[0]
		allTabs = []
		for utd in dataMap.iterUtDates():
			print 'calculating zero points for ',utd
			aperCatFn = '.'.join([pfx,utd,filt,'cat','fits'])
			files,frames = dataMap.getFiles('object',with_frames=True)
			if files is None:
				continue
			aperCat = fits.getdata(os.path.join(aperCatDir,aperCatFn))
			aperZps = np.zeros((len(frames),4),dtype=np.float32)
			psfZps = np.zeros_like(aperZps)
			for n,(f,i) in enumerate(zip(files,frames)):
				expTime =  dataMap.obsDb['expTime'][i]
				ii = np.where(aperCat['frameNum']==i)[0]
				if len(ii)==0:
					print 'no data for frame ',f
					continue
				xCat = fits.open(dataMap('cat')(f))
				for ccd in range(1,5):
					# first for the aperture photometry
					c = np.where(aperCat['ccdNum'][ii]==ccd)[0]
					mask = ( (aperCat['counts'][ii[c],aperNum]<=0) |
					         ~is_mag[aperCat['idx'][ii[c]]] )
					counts = np.ma.masked_array(
					            aperCat['counts'][ii[c],aperNum],mask=mask)
					aperMags = -2.5*np.ma.log10(counts/expTime)
					snr = counts / aperCat['countsErr'][ii[c],aperNum]
					refMags = sdss[filt][aperCat['idx'][ii[c]]]
					dMag = sigma_clip(refMags - aperMags)
					zp = np.ma.average(dMag,weights=snr**2)
					aperZps[n,ccd-1] = zp
					# then for the sextractor PSF mags
					m1,m2,s = srcor(xCat[ccd].data['ALPHA_J2000'],
					                xCat[ccd].data['DELTA_J2000'],
					                sdss['ra'][ref_ii],sdss['dec'][ref_ii],2.5)
					if len(m1)==0:
						print f,' has no catalog matches!!!'
						continue
					refMags = sdss[filt][ref_ii[m2]]
					psfMags = xCat[ccd].data['MAG_PSF'][m1]
					dMag = sigma_clip(refMags - psfMags)
					zp = np.ma.average(dMag)#,weights=snr**2)
					# have to convert from the sextractor zeropoint
					zp += 25.0 - 2.5*np.log10(expTime)
					psfZps[n,ccd-1] = zp
			tab = Table([np.repeat(utd,len(frames)),frames,aperZps,psfZps],
			            names=('utDate','frameNum','aperZp','psfZp'),
			            dtype=('S8','i4','f4','f4'))
			allTabs.append(tab)
		tab = vstack(allTabs)
		tab.write('zeropoints_%s.fits'%filt)

if __name__=='__main__':
	parser = bokrmpipe.init_file_args()
	parser.add_argument('-z','--zeropoint',action='store_true',
	                help='do zero point calculation')
	# XXX for now
	parser.add_argument('-n','--newfiles',action='store_true',
	                help='process to new files (not in-place)')
	parser.add_argument('--darkskyframes',action='store_true',
	                help='load only the dark sky frames')
	parser.add_argument('--tmpdirin',action='store_true',
	                help='read files from temporary directory')
	parser.add_argument('--tmpdirout',action='store_true',
	                help='write files to temporary directory')
	args = parser.parse_args()
	dataMap = bokrmpipe.init_data_map(args)
	if args.zeropoint:
		calc_zero_points(dataMap)

