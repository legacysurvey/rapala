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

def _get_ps1_bokmag(band,ps1m,j):
	coeffsf = os.path.join(os.environ['BOKPIPE'],'..','survey',
	                       'config','bok2ps1_%s_gicoeff.dat'%band)
	ps1colorterms = np.loadtxt(coeffsf)
	ps1mag = ps1m['MEDIAN'][:,j]
	ps1gi = np.diff(ps1m['MEDIAN'][:,[2,0]],axis=1).squeeze()
	return ps1mag + np.polyval(ps1colorterms,ps1gi)

def _residual_map(inp,ccdnum=1,nbin=8,version=None):
	filter,zp,ps1mf = inp
	print ps1mf
	try:
		nbinx,nbiny = nbin
	except:
		nbinx = nbiny = nbin
	shape = (4032//nbiny,4096//nbinx)
	if version.startswith('noao'):
		shape = shape[::-1]
	astrom_resid = np.zeros((2,)+shape,dtype=np.float32)
	astrom_nobj = np.zeros(shape,dtype=np.int32)
	photom_resid = np.zeros(shape,dtype=np.float32)
	photom_nobj = np.zeros(shape,dtype=np.int32)
	try:
		ps1m = fits.getdata(ps1mf)
	except IOError:
		return (astrom_resid,astrom_nobj,[0.0,0.0],
		        photom_resid,photom_nobj,0.0)
	ps1m = ps1m[ps1m['ccdNum']==ccdnum]
	# astrometric residuals
	cosdec = np.cos(np.radians(ps1m['DEC']))
	dra = 3600*(ps1m['RA'] - ps1m['ALPHA_J2000'])*cosdec
	dde = 3600*(ps1m['DEC'] - ps1m['DELTA_J2000'])
	sep = np.sqrt(dra**2 + dde**2)
	g = ~(sigma_clip(sep,sigma=2.5,iters=2).mask)
	xi = (ps1m['X_IMAGE']/nbinx).astype(np.int32)
	yi = (ps1m['Y_IMAGE']/nbiny).astype(np.int32)
	g &= (yi<shape[0]) & (xi<shape[1]) # remove edge sources
	np.add.at(astrom_resid[0], (yi[g],xi[g]), dra[g])
	np.add.at(astrom_resid[1], (yi[g],xi[g]), dde[g])
	np.add.at(astrom_nobj,     (yi[g],xi[g]), 1)
	astromrms = [ dra[g].std(), dde[g].std() ]
	# photometric residuals
	j = 'ugriz'.find(filter)
	flux = np.ma.array(ps1m['FLUX_APER'],mask=ps1m['FLUX_APER']<=0)
	bokmag = zp - 2.5*np.ma.log10(flux)
	ps1bokmag = _get_ps1_bokmag(filter,ps1m,j)
	dmag = bokmag - ps1bokmag
	g = ( ~(sigma_clip(dmag,sigma=2.5,iters=2).mask) &
	       (ps1m['NMAG_OK'][:,j]>2) &
	       (ps1m['MEDIAN'][:,j]>16) &
	       (ps1m['MEDIAN'][:,j]<20) &
	       (yi<shape[0]) & (xi<shape[1]) )
	if True:
		dmag -= np.ma.median(dmag[g]) # correct the zeropoint
	np.add.at(photom_resid, (yi[g],xi[g]), dmag[g])
	np.add.at(photom_nobj,  (yi[g],xi[g]), 1)
	photrms = dmag[g].std()
	return (astrom_resid,astrom_nobj,astromrms,
	        photom_resid,photom_nobj,photrms)

def merge_residual_maps(ps1matches,ccdnum,nbin,nproc,version):
	pool = multiprocessing.Pool(nproc)
	_proc = partial(_residual_map,ccdnum=ccdnum,nbin=nbin,version=version)
	res = pool.map(_proc,ps1matches)
	#res = map(_proc,ps1matches)
	astrom_resid = np.sum( [ resit[0] for resit in res ], axis=0 )
	astrom_nobj  = np.sum( [ resit[1] for resit in res ], axis=0 )
	astrom_rms  = np.array([ resit[2] for resit in res ] )
	astrom_resid = np.ma.array(astrom_resid,
	                           mask=np.tile(astrom_nobj==0,(2,1,1)))
	astrom_resid = np.ma.divide(astrom_resid,astrom_nobj)
	photom_resid = np.sum( [ resit[3] for resit in res ], axis=0 )
	photom_nobj  = np.sum( [ resit[4] for resit in res ], axis=0 )
	photom_rms  = np.array([ resit[5] for resit in res ] )
	photom_resid = np.ma.array(photom_resid,mask=astrom_nobj==0)
	photom_resid = np.ma.divide(photom_resid,photom_nobj)
	return (astrom_resid.filled(np.nan),astrom_nobj,astrom_rms,
	        photom_resid.filled(np.nan),photom_nobj,photom_rms)

def make_residual_maps(ccdsFile,outdir,nbin,nproc,byutd=False,version='naoc',
                       files_only=False):
	files = []
	ccds = Table.read(ccdsFile)
	ccds = ccds[ccds['cali_ref']=='PS1'] # only images with calibration
	if byutd:
		ccds = ccds.group_by(['filter','ccdnum','date_obs'])
	else:
		ccds = ccds.group_by(['filter','ccdnum'])
	rmsTab = []
	aprefix = 'astrom_resid_%s' % version
	pprefix = 'photom_resid_%s' % version
	for k,g in zip(ccds.groups.keys,ccds.groups):
		print 'processing ',tuple(k),len(g)
		if byutd:
			filt,ccdnum,utd = k
			utdstr = '_'+utd
		else:
			filt,ccdnum = k
			utdstr = ''
		aoutfn = '%s_%s%s_ccd%d.fits' % (aprefix,filt,utdstr,ccdnum)
		poutfn = '%s_%s%s_ccd%d.fits' % (pprefix,filt,utdstr,ccdnum)
		if ccdnum==1:
			files.append((aprefix,pprefix,utdstr))
		if files_only:
			continue
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
		zpt = g['ccdzpt'+'abcd'[ccdnum-1]] #- 2.5*np.log10(g['arawgain'])
		dat = [ (b,zp,f) for b,zp,f in zip(g['filter'],zpt,dat) ]
		resmaps = merge_residual_maps(dat,ccdnum,nbin,nproc,version)
		aresid_im,aresid_n,aresid_rms = resmaps[:3]
		presid_im,presid_n,presid_rms = resmaps[3:]
		fits.writeto(aoutfn,aresid_im,clobber=True)
		fits.writeto(aoutfn.replace('.fits','_nstar.fits'),
		             aresid_n,clobber=True)
		t = g['expnum','ccdnum','filter'].copy()
		t['raRms'] = aresid_rms[:,0]
		t['decRms'] = aresid_rms[:,1]
		t['photRms'] = presid_rms
		rmsTab.append(t)
		fits.writeto(poutfn,presid_im,clobber=True)
		fits.writeto(poutfn.replace('.fits','_nstar.fits'),
		             presid_n,clobber=True)
	if len(rmsTab) > 0:
		vstack(rmsTab).write(version+'_rms.fits',overwrite=True)
	return files

def make_plots(files,version='naoc',**kwargs):
	vmin = kwargs.get('vmin',-0.15)
	vmax = kwargs.get('vmax', 0.15)
	for aprefix,pprefix,utdstr in files:
		for pfx in [aprefix,pprefix]:
			prefix = '%s%s' % (pfx,utdstr)
			if pfx==aprefix:
				imlist = ['ra','dec']
			else:
				imlist = ['phot']
			with PdfPages(prefix+'.pdf') as pdf:
				for b in 'gr':
					for i,s in enumerate(imlist):
						fig = plt.figure(figsize=(7.25,8))
						fig.subplots_adjust(0.02,0.08,0.98,0.94,0.03,0.01)
						cax = fig.add_axes([0.1,0.03,0.8,0.03])
						for ccdnum in [1,3,2,4]:
							imf = '%s_%s%s_ccd%d.fits' % (pfx,b,utdstr,ccdnum)
							try:
								im = fits.getdata(imf)
								if pfx==aprefix: im = im[i]
							except IOError:
								continue
							if version=='noaoV0':
								im = np.flipud(im.transpose())
							if False: #'photom' in prefix:
								_vmin,_vmax = -1,1
							else:
								_vmin,_vmax = vmin,vmax
							ax = fig.add_subplot(2,2,ccdnum)
							_im = ax.imshow(im,vmin=_vmin,vmax=_vmax,
							                cmap=plt.cm.coolwarm,
							                origin='lower',
							                interpolation='nearest')
							ax.xaxis.set_visible(False)
							ax.yaxis.set_visible(False)
							if ccdnum==1:
								cb = fig.colorbar(_im,cax,
								                  orientation='horizontal')
								cb.ax.tick_params(labelsize=9)
						fig.text(0.5,0.98,'%s band %s residuals' % (b,s),
						         ha='center',va='top',size=14)
						pdf.savefig()
						plt.close()

def line_plots(version='naoc'):
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
	print 'binning is ',nbin
	files = make_residual_maps(args.input,args.outputdir,
	                           nbin,args.processes,
	                           version=args.version,byutd=args.utd,
	                           files_only=args.plotonly)
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
		make_plots(files,version=args.version,**kwargs)

