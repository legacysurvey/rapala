#!/usr/bin/env python

import os
import numpy as np
import fitsio

import matplotlib.pyplot as plt

from . import bokutil
from . import bokproc
from . import bokastrom

def check_gain_bal(fileName,badPixMaskFile=None,
                   showMean=True,showMode=True,**kwargs):
	maskFits = None
	if badPixMaskFile is not None:
		maskFits = fitsio.FITS(badPixMaskFile)
	gainBalance = bokproc.BokCalcGainBalanceFactors(save_arrays=True,
                                                 mask_map=lambda f: maskFits,
	                                                **kwargs)
	gainBalance.process_files([fileName])
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(16):
		ax = plt.subplot(4,4,i+1)
		extn = 'IM%d'%(i+1)
		pix = gainBalance.arrays[i].flatten()
		modalVal = 3*np.ma.median(pix)-2*pix.mean()
		plt.hist(pix,100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(gainBalance.skyVals[i],color='r',lw=1.2)
		if showMean:
			plt.axvline(pix.mean(),color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,
		         np.ma.median(pix),pix.mean(),
		         pix.std(),pix.min(),pix.max())

def check_sky_level(fileName,maskFile=None,statreg='ccd_central_quadrant',
                    showMean=True,showMode=True,**kwargs):
	fits = fitsio.FITS(fileName)
	maskFits = None
	if maskFile is not None:
		maskFits = fitsio.FITS(maskFile)
	statsPix = bokutil.stats_region(statreg)
	plt.figure(figsize=(9,9))
	plt.subplots_adjust(0.07,0.05,0.97,0.95,0.2,0.2)
	for i in range(4):
		ax = plt.subplot(2,2,i+1)
		extn = 'CCD%d'%(i+1)
		pix = fits[extn].read()[statsPix]
		if maskFits is not None:
			pix = np.ma.masked_array(pix,
			                         mask=maskFits[extn].read()[statsPix]>0)
		medVal,rms,pix = bokutil.array_stats(pix,method='median',
		                                     clip=True,rms=True,
		                                     retArray=True,**kwargs)
		meanVal = pix.mean()
		modalVal = 3*medVal-2*meanVal
		plt.hist(pix.flatten(),100,(pix.min(),pix.max()),
		         color='0.5',edgecolor='none',alpha=0.7,label=extn)
		plt.axvline(medVal,color='r',lw=1.2)
		if showMean:
			plt.axvline(meanVal,color='g',lw=1.2,ls='-')
		if showMode:
			plt.axvline(modalVal,color='purple',lw=1.2,ls='-.')
		ax.yaxis.set_visible(False)
		print '%4s %8d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f' % \
		        (extn,(~pix.mask).sum(),modalVal,medVal,meanVal,
		         pix.std(),pix.min(),pix.max())


_scamp_diag_head = '''<html>
<head>
 <style type="text/css">
  body {font-family: monospace}
  body {font-size: xx-small}
  td {padding-right: 10px}
 </style>
<body><table border=1>
'''
_scamp_diag_foot = '''
</table>
</body></html>
'''

def run_scamp_diag(imageFiles,**kwargs):
	tabf = open(os.path.join('scamp_diag.html'),'w')
	tabf.write(_scamp_diag_head)
	tabkeys = ['frame']+['ccd%drms'%n for n in range(1,5)]
	rowstr = ''.join([r'<th>%s</th>' % k for k in tabkeys])
	tabf.write(r'<tr>'+rowstr+r'</tr>'+'\n')
	for imFile in imageFiles:
		aheadfn = imFile.replace('.fits','.ahead')
		if not os.path.exists(aheadfn):
			print imFile,' missing'
			continue
		hdrs = bokastrom.read_headers(aheadfn)
		rowstr = r'<td>%s</td>' % os.path.basename(imFile)
		for ccdNum,hdr in enumerate(hdrs,start=1):
			rms = 3600*np.sqrt(hdr['ASTRRMS1']**2 + hdr['ASTRRMS1']**2)
			if rms < 0:
				rowstr += r'<td bgcolor=gray>%.2f</td>' % rms
				#nmissing += 1
			elif rms > 0.4:
				rowstr += r'<td bgcolor=red>%.2f</td>' % rms
				#nbad += 1
			elif rms > 0.25:
				rowstr += r'<td bgcolor=yellow>%.2f</td>' % rms
				#nmarginal += 1
			else:
				rowstr += r'<td>%.2f</td>' % rms
				#ngood += 1
		tabf.write(r'<tr>'+rowstr+r'</tr>'+'\n')
	tabf.write(_scamp_diag_foot)
	tabf.close()

