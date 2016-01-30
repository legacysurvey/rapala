#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.stats import sigma_clip

import basslog

ampOrder = [ 4,  3,  2,  1,  8,  7,  6,  5,  9, 10, 11, 12, 13, 14, 15, 16 ]

nominal_gain =  np.array(
  [ 1.3, 1.3, 1.3, 1.3, 
    1.5, 1.5, 1.3, 1.5, 
    1.4, 1.4, 1.4, 1.3, 
    1.4, 1.3, 1.4, 1.4
   ] )

nominal_rdnoise = np.array(
  [ 6.6, 6.7, 6.9, 6.3, 
    7.8, 8.2, 6.7, 6.9, 
    9.9, 6.3, 8.0, 8.6, 
    9.5, 6.0, 6.1, 6.8
  ])

gainV = np.array(
  [ 1.24556017,  1.29317832,  1.31759822,  1.28293753,  
    1.44988859, 1.52633166,  1.42589855,  1.51268101,  
    1.33969975,  1.39347458, 1.3766073 ,  1.39406121,  
    1.42733335,  1.38764536,  1.79094434, 1.45403028
  ] )

rdnoiseV = np.array(
  [  7.80125856,   8.1394453 ,   7.86576653,   7.17073345,
     9.2713232 ,  10.02771378,   9.96961975,   8.54589653,
    13.21605778,  12.06299305,  10.2200613 ,  10.66010189,
    10.41293907,   7.91155386,   8.79114914,   8.53691864 
  ] )

rm_ampmap = [ 3,1,4,2, 7,5,8,6, 10,12,9,11, 14,16,13,15 ]

pts = {'marker':'.','s':7,'c':'0.35','edgecolor':'none'}

def e_var(Nsky,j):
	return (Nsky/gainV[j]) + (rdnoiseV[j]/gainV[j])**2

def e_var_nominal(Nsky,j):
	return (Nsky/nominal_gain[j]) + (nominal_rdnoise[j]/nominal_gain[j])**2

def check_gain_from_sky(survey):
	_stats = fits.getdata('tile_stats_%s.fits' % survey,2)
	# argh, should have saved shape correctly
	stats = {}
	for k in _stats.names:
		stats[k] = _stats[k].reshape(-1,16)
	skyk = 'skyMedian'
	vark = 'skyVar'
	logSky_ref = np.log10(stats[skyk][:,0])
	gratio = gainV / gainV[0]
	outfn = 'bok90prime_gaincheck_%s' % survey
	pdf = PdfPages(outfn+'.pdf')
	for j in range(1,16):
		if j%8==1:
			fig = plt.figure(figsize=(8.0,10))
			plt.subplots_adjust(0.08,0.08,0.92,0.92,0.3,0.35)
			pnum = 1
		# just to be safe
		assert np.all(stats['ampNum'][:,j]==ampOrder[j])
		#
		ax1 = plt.subplot(8,3,pnum)
		pnum += 1 
		v = stats[skyk][:,j]/stats[skyk][:,0] 
		plt.scatter(logSky_ref,v,**pts)
		plt.axhline(np.median(v),c='r',lw=0.5)
		plt.axhline(gratio[j]**-1,c='g',lw=0.5,ls='--')
		plt.xlim(2.2,4.2)
		plt.ylim(np.median(v)-0.2,np.median(v)+0.2)
		#
		ax2 = plt.subplot(8,3,pnum)
		pnum += 1 
		v = stats[vark][:,j]/stats[vark][:,0]
		plt.scatter(logSky_ref,v,**pts)
		plt.axhline(np.median(v),c='r',lw=0.5)
		cts = np.logspace(2.0,4.5,100)
		vratio = e_var(cts,j) / e_var(cts,0)
		plt.plot(np.log10(cts),vratio,c='g',lw=0.5,ls='--')
		plt.xlim(2.2,4.2)
		plt.ylim(np.median(v)-0.4,np.median(v)+0.4)
		#
		ax3 = plt.subplot(8,3,pnum)
		pnum += 1 
		v = stats[skyk][:,j]/stats[skyk][:,0]
		plt.plot(v,c=pts['c'])
		plt.axhline(np.median(v),c='r',lw=0.5)
		plt.axhline(gratio[j]**-1,c='g',lw=0.5,ls='--')
		plt.xlim(0,len(v))
		plt.ylim(np.median(v)-0.05,np.median(v)+0.05)
		#
		for ax in [ax1,ax2,ax3]:
			for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(7)
		ax1.text(0.9,0.9,'#%d/#%d' % (stats['ampNum'][0,j],ampOrder[0]),
		         ha='right',va='top',size=10,transform=ax1.transAxes)
		if j%8==1:
			ax1.set_title('sky ratio vs. log(ADU)',size=10)
			ax2.set_title('var ratio vs. log(ADU)',size=10)
			ax3.set_title('sky ratio vs. time',size=10)
		if j==8:
			ax1.set_xlabel('log(sky,ADU)',size=8)
			ax2.set_xlabel('log(sky,ADU)',size=8)
			ax3.set_xlabel('image num.',size=8)
			pdf.savefig(fig,orientation='landscape')
			plt.close()
	pdf.savefig(fig)
	plt.close()
	pdf.close()

def check_variance_from_sky(survey):
	_stats = fits.getdata('tile_stats_%s.fits' % survey,2)
	# argh, should have saved shape correctly
	stats = {}
	for k in _stats.names:
		stats[k] = _stats[k].reshape(-1,16)
	skyk = 'skyMedian'
	vark = 'skyVar'
	outfn = 'bok90prime_varcheck_%s' % survey
	pdf = PdfPages(outfn+'.pdf')
	for j in range(16):
		if j==0:
			fig = plt.figure(figsize=(8.0,10))
			plt.subplots_adjust(0.08,0.08,0.92,0.92,0.3,0.35)
			pnum = 1
		#
		ax = plt.subplot(8,2,pnum)
		v = stats[skyk][:,j]/stats[vark][:,j]
		logSky = np.log10(stats[skyk][:,j])
		plt.scatter(logSky,v,**pts)
		cts = np.logspace(2.0,4.5,100)
		plt.plot(np.log10(cts),cts/e_var(cts,j),c='r',lw=0.5)
		plt.plot(np.log10(cts),cts/e_var_nominal(cts,j),c='g',lw=0.5)
		plt.scatter(4.0,gainV[j],marker='*',s=100,
		            edgecolor='orange',facecolor='none')
		plt.xlim(2.2,4.2)
		plt.ylim(0.75,1.65)
		for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(7)
		ax.text(0.03,0.9,'#%d'%ampOrder[j],size=9,va='top',
		        transform=ax.transAxes)
		pnum += 1 
	fig.text(0.5,0.93,r'$N_{sky}(ADU)/Var(ADU)$ vs. $log(N_{sky}(ADU)$)',
	         ha='center',va='bottom',size=11)
	pdf.savefig(fig)
	plt.close()
	pdf.close()

def check_gain_rn_from_cals(which='gain',survey='bass'):
	outfn = 'bok90prime_%s_cal%s' % (survey,which)
	fn = 'bok90_%s_char.fits' % survey
	f1 = fits.getdata(fn,1)
	f2 = fits.getdata(fn,2)
	medVal = np.ma.median(sigma_clip(f1[which],
	                                 axis=0,sig=2.2,iters=2),axis=0)
	pdf = PdfPages(outfn+'.pdf')
	for j in range(16):
		if j%8==0:
			fig = plt.figure(figsize=(8.0,10))
			plt.subplots_adjust(0.08,0.06,0.94,0.92,0.3,0.35)
			pnum = 1
			units = {'gain':'e- per ADU','rdnoise':'e-'}
			fig.text(0.03,0.5,'%s [%s]' % (which,units[which]),
			         ha='left',va='center',rotation='vertical',
			         size=12)
		#
		ax = plt.subplot(8,1,pnum)
		plt.plot(f1[which][:,j],c=pts['c'],drawstyle='steps-pre')
		plt.axhline(medVal[j],c='r',lw=1.3,alpha=0.7)
		if which=='gain':
			plt.axhline(nominal_gain[j],c='g',lw=1.3,ls='--',alpha=0.7)
		else:
			plt.axhline(nominal_rdnoise[j],c='g',lw=1.3,ls='--',alpha=0.7)
		for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(7)
		ax.text(0.03,0.9,'#%d'%ampOrder[j],size=9,va='top',
		        transform=ax.transAxes)
		plt.xlim(-5,f1[which].shape[0])
		ax.text(0.7,0.9,'%.2f'%medVal[j],va='top',color='r',size=11,
		        transform=ax.transAxes)
		ymax = f1[which][:,j].max()
		if which == 'gain':
			ymax += 0.3
			plt.ylim(1.15,ymax)
		else:
			ymax += 1.0
			plt.ylim(medVal[j]-3,ymax)
		if True:
			lastut = 0
			for _i,ut in enumerate(f2['utDate']):
				if int(ut) > lastut:# + 10:
					if j%8==0:
						dy = 0.05 if which=='gain' else 0.25
						ax.text(_i,ymax+dy,ut,size=7,
						        rotation='vertical',va='bottom',ha='center')
					ax.axvline(_i,lw=1.2,color='c',ls=':')
					lastut = int(ut)
		pnum += 1 
		if j%8==7:
			pdf.savefig(fig)
			plt.close()
	pdf.close()

def sdssrm_relative_photometry():
	import os
	import bokcat,boklog
	import fitsio
	nX,nY = 4096,4032
	nX2 = nX//2
	nY2 = nY//2
	nCCD = 4
	photdir = os.path.join(os.environ['BOK90PRIMEOUTDIR'],'catalogs_v2')
	#photdir = os.path.join('.','catalogs_v2')
	filt = 'g'
	aperNum = 3
	mag0 = 25.0
	SNRcut = 10.0
	#utd2run = bok_run_index()
	nightlyLogs = boklog.load_Bok_logs()
	#frameList = build_frame_list(filt,nightlyLogs)
	refcat = bokcat.load_targets('CFHTstars')
	catpfx = 'cfhtbright'
	#magbins = np.arange(17.,20.1,0.2)
	#mi_fix = [7,8,9]
	#mstride = 5
	magbins = np.arange(17.,20.1,1.0)
	mi_fix = [1]
	mstride = 1
#	dm_stats = np.zeros((10000,len(magbins[::mstride]),16,4),dtype=np.float32)
	dm_stats = []
	skymeas_info = fits.getdata('tile_stats_sdssrm.fits',1)
	skymeas_data = fits.getdata('tile_stats_sdssrm.fits',2)
	skyval = skymeas_data['skyMean'].reshape(-1,16).mean(axis=1)
	#refMag = refcat['g']
	#refCatMag = refcat['psfMag'][:,1]
	n = 0
	for night,utd in enumerate(sorted(nightlyLogs.keys())):
		if utd=='20131222':
			continue
		try:
			catfn = '.'.join([catpfx,utd,filt,'cat','fits'])
			_fits = fitsio.FITS(os.path.join(photdir,catfn))
		except ValueError:
			continue
		print catfn
		data = _fits[1].read()
		# map the individual entries back to their reference source
		refid = np.zeros(data.size,dtype=np.int32)
		for starNum,i1,i2 in _fits[2]['TINDEX','i1','i2'][:]:
			refid[i1:i2] = starNum
		refMag = refcat['psfMag'][refid,1]
		ref_gmr = refcat['psfMag'][refid,1] - refcat['psfMag'][refid,2]
		#
		good = ( (data['flags'][:,aperNum] == 0) &
		         (data['aperCounts'][:,aperNum] > 0) & 
		         (data['aperCountsErr'][:,aperNum] <
		           (1/SNRcut)*data['aperCounts'][:,aperNum]) )
		good = np.where(good)[0]
		ccdNums = data['ccdNum'][good]
		ampNums = (data['x'][good]//nX2).astype(np.int) + \
		          2*(data['y'][good]//nY2).astype(np.int) + \
		          4*(ccdNums-1)
		counts = data['aperCounts'][good,aperNum]
		mags = mag0 - 2.5*np.log10(counts/150.)
		errs = 1.0857*data['aperCountsErr'][good,aperNum]/counts
		refMag = refMag[good]
		ref_gmr = ref_gmr[good]
		frames = np.unique(data['frameNum'][good])
		for frame in frames:
			# index into sky measurements
			fn = nightlyLogs[utd]['fileName'][frame]
			try:
				kk = np.where((skymeas_info['utDate']==utd) &
				              (skymeas_info['fileName']==fn))[0][0]
			except:
				print utd,fn,' missing'
				continue
			for amp in range(16):
				g = np.where((data['frameNum'][good]==frame) &
				             (ampNums==amp))[0]
				if len(g)==0:
					print 'frame ',frame,amp,' not enough good stars'
					continue
				mfix = np.where((refMag[g]>18.5)&(refMag[g]<19.3))[0]
				if len(mfix)<5:
#					print 'frame ',frame,amp,' not enough zero point stars',
#					print len(mfix)
					continue
				dm = mags[g[mfix]] - refMag[g[mfix]]
				mean_dm_fix = sigma_clip(dm,sig=2.5,iters=2).mean()
				ii = np.where(refMag[g]<21.1)[0]
				dm = mags[g[ii]] - refMag[g[ii]] - mean_dm_fix
#				adat = (np.repeat(amp,len(ii)),
#				        np.repeat(skyval[kk],len(ii)),refMag[g[ii]],dm)
#				a = np.array(adat,
#	                  dtype=[('ampNum','i4'),('skyADU','f4'),
#	                         ('refMag','f4'),('dmag','f4')])
#				dm_stats.append(a)
				dm_stats.append((np.repeat(amp,len(ii)),
				                 np.repeat(skyval[kk],len(ii)),counts[g[ii]],
				                 refMag[g[ii]],ref_gmr[g[ii]],dm))
				'''
				mi = np.digitize(refMag[g],magbins) - 1
				mfix = np.where(np.in1d(mi,mi_fix))[0]
				if len(mfix)<5:
					print 'frame ',frame,amp,' not enough zero point stars'
					continue
				dm = mags[g[mfix]] - refMag[g[mfix]]
				mean_dm_fix = sigma_clip(dm,sig=2.5,iters=2).mean()
				for j,mbin in enumerate(range(0,len(magbins),mstride)):
					ii = np.where(mi==mbin)[0]
					if len(ii)<=3:
						print 'frame ',frame,amp,mbin,' not enough stars'
						continue
					dm = mags[g[ii]] - refMag[g[ii]] - mean_dm_fix
					dm = sigma_clip(dm,sig=2.5,iters=2)
					print frame,amp,mbin,len(dm),dm.mean(),dm.std()
					dm_stats[n,j,amp,:] = ( mbin, skyval[kk], 
					                        dm.mean(), dm.std() )
				'''
			n += 1
#	dm_stats = dm_stats[:n]
#	np.save('sdssrm_relphot.npy',dm_stats)
	arr = np.hstack(dm_stats)
	dm_stats = np.core.records.fromarrays(arr,
	                  dtype=[('ampNum','i4'),('skyADU','f4'),('counts','f4'),
	                         ('refMag','f4'),('ref_gmr','f4'),('dmag','f4')])
#	dm_stats = np.concatenate(dm_stats)
	fitsio.write('sdssrm_relphot.fits',dm_stats,clobber=True)

# coefficients for g_Bok - g_CFHT = a (g-r)_CFHT + b
bok2cfht = np.array([ 0.11022684, -0.07223208])

def plot_relphot(colorfix=True):
	from astrotools.idmstuff import binmean
	magbins = np.arange(17.,21.1,1.0)
	p = fits.getdata('sdssrm_relphot.fits')
	mi = np.digitize(p['refMag'],magbins) - 1
	logsky = np.log10(p['skyADU'])
	skybins = np.arange(2.5,4.51,0.20)
	dmbins = np.arange(-0.15,0.151,0.01)
	xb = skybins[:-1]+np.diff(skybins)/2
	yb = dmbins[:-1]+np.diff(dmbins)/2
	good = np.where((p['dmag']>-0.2)&(p['dmag']<0.2))[0]
	if colorfix:
		dmag = p['dmag'] - np.polyval(bok2cfht,p['ref_gmr'])
		outfn = 'sdssrm_relphot'
	else:
		dmag = p['dmag']
		outfn = 'sdssrm_relphot_nocolorfix'
	pdf = PdfPages(outfn+'.pdf')
	for amp in range(16):
		if amp%8==0:
			fig = plt.figure(figsize=(8.0,10))
			plt.subplots_adjust(0.08,0.08,0.92,0.92,0.3,0.35)
			pnum = 1
		#
		amp_ii = good[np.where(p['ampNum'][good]==amp)[0]]
		for j in range(len(magbins)-1):
			mag_ii = np.where(mi[amp_ii]==j)[0]
			n,_,_ = np.histogram2d(logsky[amp_ii[mag_ii]],
			                       dmag[amp_ii[mag_ii]],
			                       [skybins,dmbins])
			b,mdm,sdm = binmean(logsky[amp_ii[mag_ii]],
			                    dmag[amp_ii[mag_ii]],
			                    skybins,std=True,clip=True)
			#plt.figure()
			ax = plt.subplot(8,4,pnum)
#			plt.hexbin(logsky[amp_ii[mag_ii]],dmag[amp_ii[mag_ii]],
#			           gridsize=40,cmap=plt.cm.Blues)
			plt.contour(xb,yb,n.transpose(),colors='gray')
			plt.errorbar(b,mdm,sdm,c='b')
#			plt.imshow(n.transpose(),interpolation='nearest')
			plt.axhline(0,c='r')
			plt.xlim(2.2,4.6)
			plt.ylim(-0.15,0.15)
			for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(7)
			if j==0:
				ax.text(0.03,0.9,'#%d'%rm_ampmap[amp],size=9,va='top',
				        transform=ax.transAxes)
			if amp%8==0:
				ax.set_title('$%.1f < g < %.1f$' % 
				             (magbins[j],magbins[j+1]),size=11)
			pnum += 1
#		return
		if amp%8==7:
			pdf.savefig(fig)
			plt.close()
	pdf.close()

def colorcheck():
	from astrotools.idmstuff import binmean
	p = fits.getdata('sdssrm_relphot.fits')
	mags = np.unique(p['refMag'])
	gmr = np.zeros_like(mags)
	cmags = np.zeros_like(mags)
	for i,m in enumerate(mags):
		if ((i+1)%100)==0: print i+1,len(mags)
		ii = np.where(p['refMag']==m)[0]
		gmr[i] = p['ref_gmr'][ii[0]]
		cmags[i] = sigma_clip(p['dmag'][ii]).mean()
	ii = np.where((cmags>-0.15)&(cmags<0.15))[0]
	cfit = np.polyfit(gmr[ii],cmags[ii],1)
	print cfit
	plt.figure(figsize=(14,5))
	plt.subplot(121)
	plt.scatter(mags,cmags,s=5,marker='.',c='0.2',edgecolor='none')
	b,yy,dyy = binmean(mags,cmags,np.arange(17,21.5,0.35),std=True,median=True)
	plt.errorbar(b,yy,dyy,color='b')
	plt.ylim(-0.5,0.5)
	plt.subplot(122)
	plt.scatter(gmr,cmags,s=5,marker='.',c='0.2',edgecolor='none')
	b,yy,dyy = binmean(gmr,cmags,np.arange(0,1.6,0.25),std=True,median=True)
	cfit2 = np.polyfit(b[1:-1],yy[1:-1],1)
	print cfit2
	plt.errorbar(b,yy,dyy,color='b')
	plt.plot(b,np.polyval(cfit,b),c='r')
	plt.plot(b,np.polyval(cfit2,b),c='orange')
	plt.xlim(0,1.7)
	plt.ylim(-0.17,0.17)

def check_bias_levels(survey='bass'):
	from itertools import cycle
	logs = basslog.load_Bok_logs('logs/')
	outfn = '%s_bias_log' % survey
	b = fits.getdata('%s_bias_log.fits'%survey)
	medBias = np.median(b['oscanMean'],axis=0)
	medRes = np.median(b['meanResidual'],axis=0)
	medRMS = np.median(b['rmsResidual'],axis=0)
	pdf = PdfPages(outfn+'.pdf')
	nampperpage = 2
	nrows = 6
	for j in range(16):
		if j%nampperpage==0:
			fig = plt.figure(figsize=(8.0,10))
			plt.subplots_adjust(0.08,0.06,0.94,0.92,0.3,0.35)
			pnum = 1
			fig.text(0.03,0.5,'%s [%s]' % ('bias','ADU'),
			         ha='left',va='center',rotation='vertical',
			         size=12)
		#
		ax1 = plt.subplot(nrows,1,pnum)
		ax1.plot(b['oscanMean'][:,j],c=pts['c'],drawstyle='steps-pre')
		ax1.axhline(medBias[j],c='r',lw=1.3,alpha=0.7)
		pnum += 1
		#
		ax2 = plt.subplot(nrows,1,pnum)
		ax2.plot(b['meanResidual'][:,j],c=pts['c'],drawstyle='steps-pre')
		ax2.axhline(medRes[j],c='r',lw=1.3,alpha=0.7)
		ax2.axhline(0,c='gray',ls='--')
		pnum += 1 
		#
		ax3 = plt.subplot(nrows,1,pnum)
		ax3.plot(b['rmsResidual'][:,j],c=pts['c'],drawstyle='steps-pre')
		ax3.axhline(medRMS[j],c='r',lw=1.3,alpha=0.7)
		pnum += 1 
		#
		ymin,ymax = medBias[j] - 250, medBias[j] + 250
		ax1.set_ylim(ymin,ymax)
		ax2.set_ylim(medRes[j]-8,medRes[j]+10)
		ax3.set_ylim(medRMS[j]-3,medRMS[j]+5)
		for ax in [ax1,ax2,ax3]:
			for tick in ax.xaxis.get_major_ticks()+ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(7)
			ax.set_xlim(-5,b.size+5)
		for ax,t in zip([ax1,ax2,ax3],['oscan','res','rms']):
			ax.text(0.03,0.9,'#%d %s' % (ampOrder[j],t),size=9,va='top',
			        transform=ax.transAxes)
		ax1.text(0.7,0.9,'%.2f'%medBias[j],va='top',color='r',size=11,
		         transform=ax1.transAxes)
		ax2.text(0.7,0.9,'%.2f'%medRes[j],va='top',color='r',size=11,
		         transform=ax2.transAxes)
		ax3.text(0.7,0.9,'%.2f'%medRMS[j],va='top',color='r',size=11,
		         transform=ax3.transAxes)
		if True:
			uts,ii = np.unique(b['utDate'],return_index=True)
			if j%nampperpage==0:
				dy = 10
				for i,ut in zip(ii,uts):
					ax1.text(i,ymax+dy,ut,size=7,
					         rotation='vertical',va='bottom',ha='center')
			ii2 = np.concatenate([ii,[b.size]])
			for i1,i2,c in zip(ii2[:-1],ii2[1:],cycle(['white','0.7'])):
				for ax in [ax1,ax2,ax3]:
					ax.axvspan(i1,i2,color=c,alpha=0.25,ec='none')
		if True:
			ii = []
			lastType = []
			for i,ut in enumerate(b['utDate']):
				k = np.where(logs[ut]['fileName'] == b['fileName'][i])[0][0]
				if k > 0 and logs[ut]['imType'][k-1] != 'zero':
					ii.append(i)
					lastType.append(logs[ut]['imType'][k-1])
			_clr = {'object':'orange','flat':'m','dark':'0.8'}
			for i,lt in zip(ii,lastType):
				for ax in [ax1,ax2,ax3]:
					ax.axvline(i,lw=1.2,color=_clr[lt],ls=':')
		if j%nampperpage==nampperpage-1:
			pdf.savefig(fig)
			plt.close()
	pdf.close()

if __name__=='__main__':
	sdssrm_relative_photometry()

