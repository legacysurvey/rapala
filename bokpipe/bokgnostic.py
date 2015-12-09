#!/usr/bin/env python

import os
from collections import defaultdict
import numpy as np
from astropy.table import Table
import fitsio

import matplotlib.pyplot as plt

from . import bokutil
from . import bokproc
from . import bokastrom

def obs_meta_data(dataMap,outFile='obsmetadata.fits'):
	files_and_frames = dataMap.getFiles(with_frames=True)
	cols = defaultdict(list)
	for f,i in zip(*files_and_frames):
		frameId,expTime = dataMap.obsDb['frameIndex','expTime'][i]
		procf = dataMap('proc2')(f)
		catf = dataMap('cat')(f)
		try:
			imFits = fitsio.FITS(procf)
		except:
			print frameId,procf,' does not exist'
			continue
		cols['frameId'].append(frameId)
		hdrs = [ imFits[extName].read_header()
		                 for extName in ['CCD1','CCD2','CCD3','CCD4'] ]
		cols['biasDN'].append([ h['OSCANMED'] for h in hdrs ])
		cols['skyElPerSec'].append(hdrs[0]['SKYVAL']/expTime)
		try:
			catFits = fitsio.FITS(catf)
			fwhm = []
			for extNum in range(1,5):
				objs = catFits[extNum].read()
				ii = np.where((objs['FLAGS']==0) &
				              (objs['CLASS_STAR'] > 0.5))[0]
				if len(ii) > 5:
					fwhm.append(np.median(objs['FWHM_IMAGE'][ii]))
				else:
					fwhm.append(-1)
		except:
			print i,catf,' does not exist'
			fwhm = [-1,-1,-1,-1]
		cols['fwhmPix'].append(fwhm)
	cols = dict(cols)
	tab = Table(cols)
	tab.write(outFile,overwrite=True)

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


html_diag_head = '''<html>
<head>
 <style type="text/css">
  body {font-family: monospace}
  body {font-size: xx-small}
  td {padding-right: 10px}
 </style>
<body><table border=1>
'''

html_diag_foot = '''
</table>
</body></html>
'''

def html_table_entry(val,status):
	clrs = {
	  'missing':'black',
	  'nominal':'white',
	  'warning':'yellow',
	  'bad':'red',
	  'weird':'grey',
	}
	clr = clrs.get('val','white')
	return r'<td bgcolor=%s>%s</td>' % (clr,val)

def run_scamp_diag(imageFiles,ncols=4,**kwargs):
	tabf = open(os.path.join('scamp_diag.html'),'w')
	tabf.write(html_diag_head)
	rowstr = ''
	for n,imFile in enumerate(imageFiles):
		aheadfn = imFile.replace('.fits','.ahead')
		rowstr += r'<td>%s</td>' % os.path.basename(imFile)
		if not os.path.exists(aheadfn):
			print imFile,' missing'
			rowstr += html_table_entry('','missing')
		else:
			hdr = bokastrom.read_headers(aheadfn)[0]
			rms = 3600*np.sqrt(hdr['ASTRRMS1']**2 + hdr['ASTRRMS1']**2)
			if rms < 0:
				status = 'weird'
			elif rms > 0.4:
				status = 'bad'
			elif rms > 0.2:
				status = 'warning'
			else:
				status = 'nominal'
			rowstr += html_table_entry('%.2f'%rms,status)
		if (n%ncols)==ncols-1:
			tabf.write(r'<tr>'+rowstr+r'</tr>'+'\n')
			rowstr = ''
	tabf.write(html_diag_foot)
	tabf.close()

