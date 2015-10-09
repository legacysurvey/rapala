#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import bokproc

def plot_gain_vals(diagfile):
	g = np.load(diagfile)#,gains=gainCorV,skys=skyV,gainCor=gainCor)
	plt.figure(figsize=(9,6))
	plt.subplots_adjust(0.07,0.04,0.97,0.97,0.25,0.05)
	for amp in range(16):
		ax = plt.subplot(4,4,amp+1)
		plt.plot(g['gains'][:,0,amp],c='b')
		plt.axhline(g['gainCor'][0,amp],c='purple',ls='--')
		plt.plot(g['gains'][:,1,amp],c='r')
		plt.axhline(g['gainCor'][1,amp],c='orange',ls='--')
		ax.xaxis.set_visible(False)
		plt.ylim(0.91,1.09)
		ax.text(0.05,0.05,'IM%d'%bokproc.ampOrder[amp],
		        size=8,transform=ax.transAxes)
		ax.text(0.25,0.05,'%.3f'%g['gainCor'][0,amp],color='blue',
		        size=8,transform=ax.transAxes)
		ax.text(0.50,0.05,'%.3f'%g['gainCor'][1,amp],color='red',
		        size=8,transform=ax.transAxes)

