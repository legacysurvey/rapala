#!/usr/bin/env python

import os,shutil
import glob
import numpy as np
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.table import Table,vstack

from bokpipe.bokgnostic import FocalPlaneMap

def _get_ps1_bokmag(band,ps1m,j):
	coeffsf = os.path.join(os.environ['BOKPIPE'],'..','survey',
	                       'config','bok2ps1_%s_gicoeff.dat'%band)
	ps1colorterms = np.loadtxt(coeffsf)
	ps1mag = ps1m['MEDIAN'][:,j]
	ps1gi = np.diff(ps1m['MEDIAN'][:,[2,0]],axis=1).squeeze()
	return ps1mag + np.polyval(ps1colorterms,ps1gi)

def _ps1match_read_data(inp,flip=False):
	targs,b,zp,ps1mf = inp
	try:
		ps1m = fits.getdata(ps1mf)
	except IOError:
		return None
	rv = {}
	if flip:
		# for the NOAO reductions
		x = ps1m['Y_IMAGE']
		y = 4032 - ps1m['X_IMAGE']
	else:
		x = ps1m['X_IMAGE']
		y = ps1m['Y_IMAGE']
	g = (x<4096) & (y<4032)
	ps1m = ps1m[g]
	rv['x'] = x[g]
	rv['y'] = y[g]
	rv['ccdNum'] = ps1m['ccdNum']
	# index fields by filter
	key = lambda k: k+'_'+b
	# astrometric residuals
	if key('dra') in targs or key('dde') in targs:
		cosdec = np.cos(np.radians(ps1m['DEC']))
		dra = 3600*(ps1m['RA'] - ps1m['ALPHA_J2000'])*cosdec
		dde = 3600*(ps1m['DEC'] - ps1m['DELTA_J2000'])
		sep = np.sqrt(dra**2 + dde**2)
		bad = sigma_clip(sep,sigma=2.5,iters=2).mask
		if 'dra_'+b in targs:
			rv[key('dra')] = np.ma.array(dra,mask=bad).filled(np.nan)
		if 'dde_'+b in targs:
			rv[key('dde')] = np.ma.array(dde,mask=bad).filled(np.nan)
	# photometric residuals
	dmagk = key('dmag')
	if dmagk in targs:
		j = 'ugriz'.find(b)
		flux = np.ma.array(ps1m['FLUX_APER'],mask=ps1m['FLUX_APER']<=0)
		bokmag = zp - 2.5*np.ma.log10(flux)
		ps1bokmag = _get_ps1_bokmag(b,ps1m,j)
		dmag = bokmag - ps1bokmag
		g = ( ~(sigma_clip(dmag,sigma=2.5,iters=2).mask) &
		       (ps1m['NMAG_OK'][:,j]>2) &
		       (ps1m['MEDIAN'][:,j]>16) &
		       (ps1m['MEDIAN'][:,j]<20) )
		dmag.mask |= ~g
		if True:
			dmag -= np.ma.median(dmag) # correct the zeropoint
		rv[dmagk] = dmag.filled(np.nan)
	for k in ['x2','y2','xy','a','b','theta']:
		if k in targs:
			rv[key(k)] = ps1m[k.upper()+'_IMAGE']
	for k1,k2 in [('ellip','ELLIPTICITY'),('elong','ELONGATION')]:
		if k1 in targs:
			rv[key(k1)] = ps1m[k2]
	return rv

def _noao_ps1match_read_data(inp):
	return _ps1match_read_data(flip=True)

def make_residual_maps(ccdsFile,outdir,nbin,nproc,byutd=False,version='naoc',
                       doplots=False,**kwargs):
	files = []
	ccds = Table.read(ccdsFile)
	ccds = ccds[ccds['cali_ref']=='PS1'] # only images with calibration
	if byutd:
		ccds = ccds.group_by(['filter','date_obs'])
	else:
		ccds = ccds.group_by(['filter'])
	targlist = ('dra','dde','dmag','x2','y2','xy',
	            'a','b','theta','ellip','elong')
	for k,g in zip(ccds.groups.keys,ccds.groups):
		print 'processing ',tuple(k),len(g)
		filt = str(k['filter'])
		if byutd:
			utdstr = str(k['utd'])
			targs = [ '_'.join([t,filt,utdstr]) for t in targlist ]
		else:
			targs = [ '_'.join([t,filt]) for t in targlist ]
		if version=='naoc':
#			dat = [ os.path.join(outdir,'p%s%s%s_%d.ps1match.fits') % 
#			                               (expstr[:4],filt,expstr[4:8],ccdnum) 
#			            for expstr in g['expnum'].astype(str) ]
			dat = [ os.path.join(outdir,
			           os.path.basename(f).replace('.fits','.ps1match.fits'))
			            for f in g['image_filename'] ]
		elif version=='noaoV0':
			dat = [ os.path.join(outdir,
			           'ksb_%s_%s_ooi_%s_v0.ps1match.fits' %
			              (d[2:].replace('-',''),u[:8].replace(':',''),f))
			            for d,u,f in g['date_obs','ut','filter'] ]
		zpt = np.choose(g['ccdnum']-1,[g['ccdzpt'+n] for n in 'abcd'])
		dat = [ (targs,filt,zp,f) for zp,f in zip(zpt,dat) ]
		bok2ps1fpMap = FocalPlaneMap(nbin,targs,_ps1match_read_data,nproc,
		                             prefix='fpmap_'+version)
		bok2ps1fpMap.ingest(dat)
		bok2ps1fpMap.write()
		if doplots:
			bok2ps1fpMap.make_plots(**kwargs)

def line_plots(version='naoc'):
	raise NotImplementedError
	if version=='naoc':
		ims = np.array([fits.getdata('astrom_resid_naoc_r_ccd%d.fits'%i) 
		                       for i in range(1,5)])
		nobj = np.array([fits.getdata('astrom_resid_naoc_r_ccd%d_nstar.fits'%i) 
		                       for i in range(1,5)])
	elif version=='noaoV0':
		ims = np.array([fits.getdata('astrom_resid_noaoV0_r_ccd%d.fits'%i) 
		                       for i in range(1,5)])
		ims = ims.swapaxes(2,3)[:,:,::-1,:]
		nobj = np.array(
		        [fits.getdata('astrom_resid_noaoV0_r_ccd%d_nstar.fits'%i) 
		                       for i in range(1,5)])
		nobj = nobj.swapaxes(2,3)[:,:,::-1,:]
	ims = np.ma.array(ims,mask=np.tile(nobj==0,(2,1,1,1)).swapaxes(0,1))
	nbiny = 4032 // ims.shape[2]
	nbinx = 4096 // ims.shape[3]
	yi = np.arange(ims.shape[2])*nbiny
	xi = np.arange(ims.shape[3])*nbinx
	ymid = ims.shape[2]//2
	xmid = ims.shape[3]//2
	d1 = sigma_clip(ims[:,0,ymid-1]-ims[:,0,ymid],axis=1).mean(axis=1)
	d2 = sigma_clip(ims[:,1,ymid-1]-ims[:,1,ymid],axis=1).mean(axis=1)
	d3 = sigma_clip(ims[:,0,:,xmid-1]-ims[:,0,:,xmid],axis=1).mean(axis=1)
	d4 = sigma_clip(ims[:,1,:,xmid-1]-ims[:,1,:,xmid],axis=1).mean(axis=1)
	for i,d in enumerate([d1,d2,d3,d4]):
		s = ['d(ra) at dec boundary',
		     'd(dec) at dec boundary',
		     'd(ra) at ra boundary',
		     'd(dec) at ra boundary'][i]
		print '%25s  ' % s,
		print ' '.join(['%6.3f'%x for x in d])
	plt.figure()
	#
	plt.subplot(221)
	#plt.fill_between(xi,ims[:,0].std(axis=1).transpose())
	lines = plt.plot(xi,ims[:,0].mean(axis=1).transpose())
	plt.xlim(-5,xi.max()+5)
	plt.ylim(-0.1,0.1)
	plt.axvline(xmid,c='0.2',ls='--')
	plt.legend(lines,['CCD%d'%i for i in range(1,5)],ncol=2,fontsize=10)
	plt.title(r'$\Delta(ra)$ along ra axis',size=11)
	#
	plt.subplot(222)
	plt.plot(yi,ims[:,1].mean(axis=2).transpose())
	plt.axvline(ymid,c='0.2',ls='--')
	plt.xlim(-5,yi.max()+5)
	plt.ylim(-0.1,0.1)
	plt.title(r'$\Delta(dec)$ along dec axis',size=11)
	#
	plt.subplot(223)
	plt.plot(xi,ims[:,1].mean(axis=1).transpose())
	plt.axvline(xmid,c='0.2',ls='--')
	plt.xlim(-5,xi.max()+5)
	plt.ylim(-0.1,0.1)
	plt.title(r'$\Delta(dec)$ along ra axis',size=11)
	#
	plt.subplot(224)
	plt.plot(yi,ims[:,0].mean(axis=2).transpose())
	plt.axvline(ymid,c='0.2',ls='--')
	plt.xlim(-5,yi.max()+5)
	plt.ylim(-0.1,0.1)
	plt.title(r'$\Delta(ra)$ along dec axis',size=11)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("input",type=str,
	                    help="ccds file")
	parser.add_argument("-o","--outputdir",type=str,default='./',
	                    help="output directory")
	parser.add_argument("-p","--processes",type=int,default=1,
	                    help="number of processes")
	parser.add_argument("--nbin",type=str,default='32',
	                    help="bin size for residual images")
	parser.add_argument("-V","--version",type=str,default='naoc',
	                    help="pipeline version")
	parser.add_argument("--plots",action="store_true",
	                    help="make plots")
	parser.add_argument("--range",type=str,
	                    help="plot range")
	parser.add_argument("--utd",action="store_true",
	                    help="by utdate")
	parser.add_argument("--plotonly",action="store_true",
	                    help="only make plot")
	args = parser.parse_args()
	nbin = [int(v) for v in args.nbin.split(',')]
	if len(nbin)==1: nbin = nbin[0]
	print 'binning is ',nbin
	if args.plots:
		kwargs = {}
		if args.range:
			try:
				vmin,vmax = [float(v) for v in args.range.split(',')]
			except:
				vmax = float(args.range)
				vmin = -vmax
			kwargs['vmin'] = vmin
			kwargs['vmax'] = vmax
	fpmap = make_residual_maps(args.input,args.outputdir,
	                           nbin,args.processes,
	                           version=args.version,byutd=args.utd,
	                           doplots=args.plots,**kwargs)

