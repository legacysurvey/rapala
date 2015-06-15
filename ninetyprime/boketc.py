#!/usr/bin/env python

from math import sqrt,pi

# extinction parameters
k_ext = {'g':0.17,'r':0.10}
A_ext = {'g':3.303,'r':2.285}
#r_half = 0.45
r_half = 0.35

# DESI parameters
mag_lim = {'g':24.0,'r':23.6}

# Bok parameters
p = 0.455
G = 1.375
bok_rdnoise = 7.3
bok_medianFWHM = {'g':1.5,'r':1.5}
bok_zpt0_am13 = {'g':25.55,'r':25.23}
bok_zpt0_am10 = {k: v-k_ext[k]*(1.3-1) for k,v in bok_zpt0_am13.items()}
bok_zpt0_am00 = {k: v-k_ext[k]*1.3 for k,v in bok_zpt0_am13.items()}
kpno_sky_lun0 = {'g':22.10,'r':21.07}

def _exp_config(band,mag=None,zpt0=None,sky=None,fwhm=None,
	            profile='exponential',verbose=False,tablezpt=False):
	if mag is None:
		mag = mag_lim[band]
	if zpt0 is None:
		if tablezpt:
			zpt0 = bok_zpt0_am13[band]
		else:
			zpt0 = bok_zpt0_am10[band]
	if sky is None:
		sky = kpno_sky_lun0[band]
	if fwhm is None:
		fwhm = bok_medianFWHM[band]
	#
	f0 = 10**(-0.4*(mag - zpt0)) * G
	#skyflux = 10**(-0.4*(sky - zpt0)) * p**2 * G
	skyflux = 6.6 # from images
	sig_sky = sqrt(skyflux)
	sig_readout = bok_rdnoise
	sig_seeing = fwhm / 2.355
	#
	if profile=='exponential':
		NEA = 4 * pi * ((sig_seeing/p)**2 + p**2/12 + (r_half/p)**2)
	elif profile=='gaussian':
		NEA = 4 * pi * (sig_seeing/p)**2
	else:
		raise ValueError
	if verbose:
		WNEG = 2.35*sqrt(NEA/(4*pi))
		print 'zero point = %.2f :: ' % zpt0,
		print 'flux=%.2f, sky=%.2f, rms=%.2f, NEA=%.2f, NEG=%.2f' % \
		            (f0,skyflux,sig_sky,NEA,WNEG*p)
	return f0,skyflux,sig_sky,NEA,sig_readout

def targalt_eqn1(band,SNR=6.,ndith=3,**kwargs):
	f0,skyflux,sig_sky,NEA,sig_readout = _exp_config(band,**kwargs)
	#
	T_exp = SNR**2 * NEA * (sig_sky/f0)**2 * \
	         ( 1./2 + sqrt( 1./4 + ndith/(SNR**2 * NEA * (sig_sky/f0)**2) *
	                                    (sig_readout/sig_sky)**2 ) )
	if kwargs.get('verbose',False):
		print 'calculated exposure time: %.2f' % T_exp
		print 'source counts: %.1f e-  %.1f ADU' % (f0*T_exp,f0*T_exp/G)
		print 'sky counts: %.1f e-  %.1f ADU' % (skyflux*T_exp,skyflux*T_exp/G)
	return T_exp

def snr_singleexposure(band,mag,texp,**kwargs):
	f0,skyflux,sig_sky,NEA,sig_readout = _exp_config(band,mag=mag,**kwargs)
	snr = ( f0*texp / 
	             sqrt( NEA * (skyflux*texp + sig_readout**2) ) )
	if kwargs.get('verbose',False):
		print 'calculated SNR: %.2f' % snr
		print 'source counts: %.1f e-  %.1f ADU' % (f0*texp,f0*texp/G)
		print 'sky counts: %.1f e-  %.1f ADU' % (skyflux*texp,skyflux*texp/G)
	return snr

def texp_onsky(band,airmass,ebv,mag=None,**kwargs):
	if mag is None:
		mag = mag_lim[band]
	mag += k_ext[band]*airmass + A_ext[band]*ebv
	return targalt_eqn1(band,mag=mag,**kwargs)

