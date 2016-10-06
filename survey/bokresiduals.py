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

def _residual_map(ps1mf,nbin=8,version=None):
	print ps1mf
	shape = (4032//nbin,4096//nbin)
	if version.startswith('noao'):
		shape = shape[::-1]
	astrom_resid = np.zeros((2,)+shape,dtype=np.float32)
	astrom_nobj = np.zeros(shape,dtype=np.int32)
	try:
		ps1m = fits.getdata(ps1mf)
	except IOError:
		rms = [ 0.0, 0.0 ]
		return (astrom_resid,astrom_nobj,rms)
	cosdec = np.cos(np.radians(ps1m['DEC']))
	dra = 3600*(ps1m['RA'] - ps1m['ALPHA_J2000'])*cosdec
	dde = 3600*(ps1m['DEC'] - ps1m['DELTA_J2000'])
	sep = np.sqrt(dra**2 + dde**2)
	g = ~(sigma_clip(sep,sigma=2.5,iters=2).mask)
	xi = (ps1m['X_IMAGE']/nbin).astype(np.int32)
	yi = (ps1m['Y_IMAGE']/nbin).astype(np.int32)
	g &= (yi<shape[0]) & (xi<shape[1]) # remove edge sources
	np.add.at(astrom_resid[0], (yi[g],xi[g]), dra[g])
	np.add.at(astrom_resid[1], (yi[g],xi[g]), dde[g])
	np.add.at(astrom_nobj, (yi[g],xi[g]), 1)
	rms = [ dra[g].std(), dde[g].std() ]
	return (astrom_resid,astrom_nobj,rms)

def merge_residual_maps(ps1matches,nbin,nproc,version):
	pool = multiprocessing.Pool(nproc)
	_proc = partial(_residual_map,nbin=nbin,version=version)
	res = pool.map(_proc,ps1matches)
	astrom_resid = np.sum( [ resid for resid,nobj,rms in res ], axis=0 )
	astrom_nobj  = np.sum( [ nobj  for resid,nobj,rms in res ], axis=0 )
	astrom_rms  = np.array( [ rms  for resid,nobj,rms in res ] )
	astrom_resid = np.ma.array(astrom_resid,
	                           mask=np.tile(astrom_nobj==0,(2,1,1)))
	astrom_resid = np.ma.divide(astrom_resid,astrom_nobj)
	return astrom_resid.filled(np.nan),astrom_nobj,astrom_rms

def make_residual_maps(ccdsFile,outdir,nbin,nproc,byutd=False,version='naoc'):
	ccds = Table.read(ccdsFile)
	ccds = ccds[ccds['cali_ref']=='PS1'] # only images with calibration
	if byutd:
		ccds = ccds.group_by(['filter','ccdnum','date_obs'])
	else:
		ccds = ccds.group_by(['filter','ccdnum'])
	rmsTab = []
	prefix = 'astrom_resid_%s' % version
	for k,g in zip(ccds.groups.keys,ccds.groups):
		print 'processing ',k,len(g)
		if byutd:
			filt,ccdnum,utd = k
			utdstr = '_'+utd
		else:
			filt,ccdnum = k
			utdstr = ''
		outfn = '%s_%s_ccd%d%s.fits'%(prefix,filt,ccdnum,utdstr)
		if os.path.exists(outfn):
			continue
		if version=='naoc':
#			files = [ os.path.join(outdir,'p%s%s%s_%d.ps1match.fits') % 
#			                               (expstr[:4],filt,expstr[4:8],ccdnum) 
#			            for expstr in g['expnum'].astype(str) ]
			files = [ os.path.join(outdir,
			           os.path.basename(f).replace('.fits','.ps1match.fits'))
			            for f in g['image_filename'] ]
		elif version=='noaoV0':
			files = [ os.path.join(outdir,
			           'ksb_%s_%s_ooi_%s_v0.ps1match.fits' %
			              (d[2:].replace('-',''),u[:8].replace(':',''),f))
			            for d,u,f in g['date_obs','ut','filter'] ]
		resid_im,resid_n,resid_rms = merge_residual_maps(files,nbin,nproc,
		                                                 version)
		fits.writeto(outfn,resid_im,clobber=True)
		fits.writeto(outfn.replace('.fits','_nstar.fits'),resid_n,clobber=True)
		t = g['expnum','ccdnum','filter'].copy()
		t['raRms'] = resid_rms[:,0]
		t['decRms'] = resid_rms[:,1]
		rmsTab.append(t)
	vstack(rmsTab).write(prefix+'rms.fits')

def make_plots(prefix='astrom_resid',**kwargs):
	vmin = kwargs.get('vmin',-0.15)
	vmax = kwargs.get('vmax',-0.15)
	with PdfPages(prefix+'.pdf') as pdf:
		for b in 'gr':
			for i,s in enumerate(['ra','dec']):
				fig = plt.figure(figsize=(7.25,8))
				fig.subplots_adjust(0.02,0.08,0.98,0.94,0.03,0.01)
				cax = fig.add_axes([0.1,0.03,0.8,0.03])
				for ccdnum in [1,3,2,4]:
					imf = '_'.join([prefix,b,'ccd%d.fits'%ccdnum])
					im = fits.getdata(imf)
					ax = fig.add_subplot(2,2,ccdnum)
					_im = ax.imshow(im[i],vmin=vmin,vmax=vmax,
					                cmap=plt.cm.coolwarm,
					                origin='lower',interpolation='nearest')
					ax.xaxis.set_visible(False)
					ax.yaxis.set_visible(False)
					if ccdnum==1:
						cb = fig.colorbar(_im,cax,orientation='horizontal')
						cb.ax.tick_params(labelsize=9)
				fig.text(0.5,0.98,'%s band %s residuals' % (b,s),
				         ha='center',va='top',size=14)
				pdf.savefig()
				plt.close()

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("input",type=str,
	                    help="ccds file")
	parser.add_argument("-o","--outputdir",type=str,default='./',
	                    help="output directory")
	parser.add_argument("-p","--processes",type=int,default=1,
	                    help="number of processes")
	parser.add_argument("--nbin",type=int,default=32,
	                    help="bin size for residual images")
	parser.add_argument("-V","--version",type=str,default='naoc',
	                    help="pipeline version")
	parser.add_argument("--plots",action="store_true",
	                    help="make plots")
	parser.add_argument("--range",type=str,
	                    help="plot range")
	args = parser.parse_args()
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
		make_plots(prefix=args.input,**kwargs)
	else:
		make_residual_maps(args.input,args.outputdir,
		                   args.nbin,args.processes,
		                   version=args.version)

