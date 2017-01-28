#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc
from astropy.stats import sigma_clip
import astropy.visualization as vis
from astropy.visualization.mpl_normalize import ImageNormalize
import fitsio

from . import bokutil

def make_fov_image(fov,pngfn=None,**kwargs):
	stretch = kwargs.get('stretch','linear')
	interval = kwargs.get('interval','zscale')
	imrange = kwargs.get('imrange')
	contrast = kwargs.get('contrast',0.25)
	ccdplotorder = ['CCD2','CCD4','CCD1','CCD3']
	if interval == 'rms':
		try:
			losig,hisig = imrange
		except:
			losig,hisig = (2.5,5.0)
	#
	cmap = kwargs.get('cmap','jet')
	cmap = plt.get_cmap(cmap)
	cmap.set_bad('w',1.0)
	w = 0.4575
	h = 0.455
	rc('text',usetex=False)
	fig = plt.figure(figsize=(6,6.5))
	cax = fig.add_axes([0.1,0.04,0.8,0.01])
	ims = [ fov[ccd]['im'] for ccd in ccdplotorder ]
	allpix = np.ma.array(ims).flatten()
	stretch = {
	  'linear':vis.LinearStretch(),
	  'histeq':vis.HistEqStretch(allpix),
	  'asinh':vis.AsinhStretch(),
	}[stretch]
	if interval=='zscale':
		iv = vis.ZScaleInterval(contrast=contrast)
		vmin,vmax = iv.get_limits(allpix)
	elif interval=='rms':
		nsample = 1000 // nbin
		background = sigma_clip(allpix[::nsample],iters=3,sigma=2.2)
		m,s = background.mean(),background.std()
		vmin,vmax = m-losig*s,m+hisig*s
	elif interval=='fixed':
		vmin,vmax = imrange
	else:
		raise ValueError
	norm = ImageNormalize(vmin=vmin,vmax=vmax,stretch=stretch)
	for n,(im,ccd) in enumerate(zip(ims,ccdplotorder)):
		if im.ndim == 3:
			im = im.mean(axis=-1)
		x = fov[ccd]['x']
		y = fov[ccd]['y']
		i = n % 2
		j = n // 2
		pos = [ 0.0225 + i*w + i*0.04, 0.05 + j*h + j*0.005, w, h ]
		ax = fig.add_axes(pos)
		_im = ax.imshow(im,origin='lower',
		                extent=[x[0,0],x[0,-1],y[0,0],y[-1,0]],
		                norm=norm,cmap=cmap,
		                interpolation=kwargs.get('interpolation','nearest'))
		if fov['coordsys']=='sky':
			ax.set_xlim(x.max(),x.min())
		else:
			ax.set_xlim(x.min(),x.max())
		ax.set_ylim(y.min(),y.max())
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
		if n == 0:
			cb = fig.colorbar(_im,cax,orientation='horizontal')
			cb.ax.tick_params(labelsize=9)
	tstr = fov.get('file','')+' '+fov.get('objname','')
	title = kwargs.get('title',tstr)
	title = title[-60:]
	fig.text(0.5,0.99,title,ha='center',va='top',size=12)
	if pngfn is not None:
		plt.savefig(pngfn)
		plt.close(fig)

def make_fov_image_fromfile(fileName,pngfn,nbin=1,coordsys='sky',**kwargs):
	fits = bokutil.BokMefImage(fileName,
	                           mask_file=kwargs.get('mask'),
	                           mask_type=kwargs.get('mask_type'),
	                           read_only=True)
	fov = fits.make_fov_image(nbin,coordsys)
	fov['file'] = fileName
	return make_fov_image(fov,pngfn,**kwargs)

