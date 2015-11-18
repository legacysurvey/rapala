#!/usr/bin/env python

import numpy as np
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc

from . import bokutil

def make_fov_image(fov,pngfn=None,**kwargs):
	maskFile = kwargs.get('mask')
	losig = kwargs.get('lo',2.5)
	hisig = kwargs.get('hi',5.0)
	cmap = kwargs.get('cmap','jet')
	print 'using colormap ',cmap
	cmap = plt.get_cmap(cmap)
	cmap.set_bad('w',1.0)
	w = 0.4575
	h = 0.455
	if maskFile is not None:
		maskFits = fitsio.FITS(maskFile)
	input_vmin = kwargs.get('vmin')
	input_vmax = kwargs.get('vmax')
	rc('text',usetex=False)
	fig = plt.figure(figsize=(6,6.5))
	cax = fig.add_axes([0.1,0.04,0.8,0.01])
	for n,ccd in enumerate(['CCD2','CCD4','CCD1','CCD3']):
		im = fov[ccd]['im']
		if maskFile is not None:
			im = np.ma.masked_array(im,maskFits[ccd][:,:].astype(bool))
		if n == 0:
			i1,i2 = 100//fov['nbin'],1500//fov['nbin']
			if input_vmin is None and input_vmax is None:
				background = sigma_clip(im[i1:i2,i1:i2],iters=3,sig=2.2)
				m,s = background.mean(),background.std()
				vmin = input_vmin if input_vmin is not None else m-losig*s
				vmax = input_vmax if input_vmax is not None else m+hisig*s
			else:
				vmin = input_vmin
				vmax = input_vmax
			norm = colors.Normalize(vmin=vmin,vmax=vmax)
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
	title = kwargs.get('title',fov.get('file','')+' '+fov.get('objname',''))
	fig.text(0.5,0.99,title,ha='center',va='top',size=12)
	if pngfn is not None:
		plt.savefig(pngfn)
		plt.close(fig)

def make_fov_image_fromfile(fileName,pngfn,nbin=1,coordsys='sky',**kwargs):
	fits = bokutil.BokMefImage(fileName,mask_file=kwargs.get('mask'),
	                           read_only=True)
	fov = fits.make_fov_image(nbin,coordsys)
	fov['file'] = fileName
	return make_fov_image(fov,pngfn,**kwargs)

