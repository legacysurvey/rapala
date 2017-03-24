#!/usr/bin/env python

import os,sys
import re
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack,join,unique
from astropy.time import Time,TimeDelta

ampNums = [ [ 4,2,3,1 ] , [ 7,5,8,6 ], [ 10,12,9,11 ], [ 13,15,14,16] ]

def get_amp_index(x,y):
	nx = 4096 // 2
	ny = 4032 // 2
	xi = (x/nx).astype(np.int32)
	yi = (y/ny).astype(np.int32)
	ampIndex = 2*yi + xi
	return ampIndex

try:
	bass_dir = os.environ['BASSDIR']
	bass_data_dir = os.environ['BASSDATA']
except:
	print 'env variables BASSDIR and BASSDATA, many things will fail'
	bass_dir = None
	bass_data_dir = None

try:
	rdxdir = os.environ['BASSRDXDIR']
except:
	try:
		rdxdir = os.path.join(os.environ['GSCRATCH'],'rmreduce')
	except:
		rdxdir = None

tiledb_file = 'bass-newtiles-indesi.fits'
obsdb_file = 'bass-newtiles-observed.fits'

# filenames get written in weird ways
def reform_filename(s):
	s1,s2 = re.match('.*\w(\d\d\d\d)[\w.]+(\d\d\d\d)',s).groups()
	return 'd'+s1+'.'+s2

def get_obsdb_filename(which,newest):
	obsdbfn = obsdb_file
	if which=='all':
		obsdbfn = obsdbfn.replace('.fits','_all.fits')
	if newest:
		obsdbfn = obsdbfn.replace('.fits','_updated.fits')
	return obsdbfn

def build_obsdb(update=False,which='good',newest=True):
	'''which is "good" [only good tiles] or "all" [all observations]
	   if newest is True include most recent observations 
	      not yet in "good" lists
	'''
	import glob,re,shutil
	from urllib2 import urlopen
	import tarfile
	# update the local observations database by downloading the master
	# database from the wiki site
	if update:
		resp = urlopen('http://batc.bao.ac.cn/BASS/lib/exe/fetch.php?media=observation:observation:database.tar.gz')
		f = resp.read()
		tarname = os.path.join(bass_dir,'database.tar.gz')
		outf = open(tarname,'wb')
		outf.write(f)
		outf.close()
		shutil.rmtree(os.path.join(bass_dir,'database'))
		tar = tarfile.open(tarname)
		tar.extractall(path=bass_dir)
		tar.close()
	# the summary files listing the tiles marked as "good"
	good_files = ['obsed-g-2015-good.txt','obsed-r-2015-good.txt',
	              'obsed-g-2016-0102-good.txt','obsed-r-2016-0102-good.txt',
	              'obsed-g-2016-03-good.txt','obsed-r-2016-03-good.txt',
	              'obsed-g-2016-04-good.txt','obsed-r-2016-04-good.txt',
	              'obsed-g-2016-05-good.txt','obsed-r-2016-05-good.txt',
	]
	obsfiles_good = [os.path.join(bass_dir,'database',f) for f in good_files]
	nov15_files = glob.glob(os.path.join(bass_dir,'database',
	                                     'obsed-[gr]-2015-11-??.txt'))
	# XXX not actually sure these have been classified but incl. here for now
	obsfiles_good += nov15_files
	# the original nightly tile lists archived in "*_old" directories
	obsfiles_old = glob.glob(os.path.join(bass_dir,'database','201?_old',
	                                      'obsed-[gr]-????-??-??.txt'))
	# the most recently observed tiles that have not been ingested into
	# the "good" lists yet
	obsfiles_new = glob.glob(os.path.join(bass_dir,'database',
	                                      'obsed-[gr]-201[67]-??-??.txt'))
	# select which observations to use
	if which=='all':
		# include all the observed tiles
		obsfiles = obsfiles_old
	elif which=='good':
		# only include "good" tiles
		obsfiles = obsfiles_good
	else:
		raise ValueError
	if newest and len(obsfiles_new)>0:
		# add in the most recent observations (if there are any)
		obsfiles = sorted(obsfiles + obsfiles_new)
	else:
		newest = False
	# convert the input databases to a single FITS table with some added
	# fields
	obsdb = []
	for obsfile in obsfiles:
		# for some reason astropy.Table barfs on reading this in directly
		# so working around it
		def idconv(s):
			try:
				return int(s)
			except:
				return -99
		arr = np.loadtxt(obsfile,dtype=[('fileName','S10'),('expTime','f4'),
		                           ('tileId','i4'),('ra','f8'),('dec','f8')],
		                 converters={0:reform_filename,2:idconv})
		if arr.size<=1:
			# for some reason len() freaks out in this case
			continue
		if '2015-good' in obsfile:
			# each line in this file is for a single CCD
			arr = arr[::4]
		print obsfile,len(arr)
		t = Table(arr)
		t['ditherId'] = t['tileId'] % 10
		t['tileId'] //= 10
		t['filter'] = os.path.basename(obsfile)[6]
		# filename is encoded with last 4 digits of JD
		t['mjd'] = 50000. + np.array([int(d[1:5]) for d in arr['fileName']],
		                             dtype=np.float32)
		obsdb.append(t)
	obsdb = vstack(obsdb)
	obsdb.sort('fileName')
	#
	print 'ingested %d observed tiles' % len(obsdb)
	obsdbfn = get_obsdb_filename(which,newest)
	obsdb.write(os.path.join(bass_dir,obsdbfn),overwrite=True)
	return

def load_tiledb():
	return fits.getdata(os.path.join(bass_dir,tiledb_file))

def load_obsdb(dbfile=obsdb_file):
	return fits.getdata(os.path.join(bass_dir,dbfile))

def obsdbs_joined(newest=False):
	goodobs = Table(load_obsdb(get_obsdb_filename('good',newest)))
	allobs = Table(load_obsdb(get_obsdb_filename('all',newest)))
	goodobs['good'] = True
	return join(allobs,goodobs['fileName','good'],
	            join_type='outer',keys='fileName')

def files2tiles(obsdb,fileNames):
	idxs = { row['fileName']:i for i,row in enumerate(obsdb) }
	return np.array([idxs.get(fn,-1) for fn in fileNames])

def region_tiles(ra1,ra2,dec1,dec2,observed=True):
	if observed:
		tiledb = load_obsdb()
		ii = np.where((tiledb['ra']>ra1) & (tiledb['ra']<ra2) &
		              (tiledb['dec']>dec1) & (tiledb['dec']<dec2))[0]
	else:
		tiledb = load_tiledb()
		ii = np.where((tiledb['TRA']>ra1) & (tiledb['TRA']<ra2) &
		              (tiledb['TDEC']>dec1) & (tiledb['TDEC']<dec2))[0]
	return tiledb[ii]

def get_coverage(obsdb,tiledb):
	tileCov = np.zeros((len(tiledb),2,3),dtype=np.int32)
	tid = np.array([int(tid) for tid in tiledb['TID']])
	for n,row in enumerate(obsdb):
		if row['tileId']>0:
			try:
				i = np.where(tid==row['tileId'])[0][0]
			except:
				print 'tile ',row['tileId'],' is not in db'
				continue
			if row['filter']=='g':
				tileCov[i,0,row['ditherId']-1] += 1
			else:
				tileCov[i,1,row['ditherId']-1] += 1
	return tileCov

def obs_summary(which='good',newest=True,tiles=None,
                mjdstart=None,mjdend=None,
                doplot=False,smallplot=False,saveplot=None,
                decalsstyle=False,byfilter=False,brightcolors=False,
                verbose=0):
	from collections import defaultdict
	tiledb = load_tiledb()
	obsdb = load_obsdb(get_obsdb_filename(which,newest))
	#
	if mjdstart is not None:
		obsdb = obsdb[obsdb['mjd']>=mjdstart]
	if mjdend is not None:
		obsdb = obsdb[obsdb['mjd']<=mjdend]
	#
	tid = np.array([int(tid) for tid in tiledb['TID']])
	nobs = np.zeros((tiledb.size,2,3),dtype=int)
	tileList = {1:defaultdict(list),2:defaultdict(list),3:defaultdict(list)}
	for n,row in enumerate(obsdb):
		if row['tileId']>0:
			try:
				i = np.where(tid==row['tileId'])[0][0]
			except:
				if verbose > 0:
					print 'tile ',row['tileId'],' is not in db'
				continue
			if tiles is not None and row['tileId'] not in tiles:
				continue
			tileList[row['ditherId']][row['tileId']].append(n)
			if row['filter']=='g':
				nobs[i,0,row['ditherId']-1] += 1
			else:
				nobs[i,1,row['ditherId']-1] += 1
	tileCov = nobs > 0
	if tiles is None:
		nTiles = float(len(tid))
		dec3270 = np.where((tiledb['TDEC']>=32)&(tiledb['TDEC']<=70))[0]
		nTiles3270 = float(len(dec3270))
	else:
		nTiles3270 = nTiles = float(len(tiles))
	print
	print '  MJDs %d to %d' % (obsdb['mjd'].min(),obsdb['mjd'].max())
	print
	print ' '*5,'g band'.center(32,'-'),'   ','r band'.center(32,'-')
	print ' '*5,
	print '%5s  %5s  %8s  %8s    ' % ('total','uniq','%compl','%(32-70)'),
	print '%5s  %5s  %8s  %8s' % ('total','uniq','%compl','%(32-70)')
	for _j in range(4):
		if _j < 3:
			print ' P%d: ' % (_j+1),
			j,_n = _j,1.0
		else:
			print 'all: ',
			j,_n = slice(None),3.0
		for i,filt in enumerate('gr'):
			print '%5d ' % (np.sum(nobs[:,i,j])),
			print '%5d ' % (np.sum(tileCov[:,i,j])),
			print '%8.1f ' % (100*np.sum(tileCov[:,i,j])/nTiles/_n),
			print '%8.1f ' % (100*np.sum(tileCov[:,i,j])/nTiles3270/_n),
			print '  ',
		print
	print
	if tiles is not None and verbose > 0:
		print 'MJDs of completed tiles (last MJD shown when repeated)'
		print '%6s    %s' % ('tile','g band'.center(23,'-')),
		print ' %s' % ('r band'.center(23,'-'))
		print '%6s    %7s %7s %7s %7s %7s %7s' % \
		        tuple(['']+['P%d'%(_p) 
		                for _b in 'gr' for _p in range(1,4)])
		_mjds = []
		for t in tiles:
			i = np.where(t==tid)[0][0]
			print '%6d   ' % t,
			for j,b in enumerate('gr'):
				for k in range(3):
					if tileCov[i,j,k]:
						kk = np.where((obsdb['tileId']==t)&
						              (obsdb['filter']==b)&
						              (obsdb['ditherId']==k+1))[0]
						_mjds.append(obsdb['mjd'][kk[-1]])
						print '%7d' % obsdb['mjd'][kk[-1]],
					else:
						print '%7s' % '---',
			print
		print
		uts = [ Time(mjd,format='mjd').iso[:10] 
		           for mjd in sorted(np.unique(_mjds)) ]
		print 'unique set of UT dates covering these tiles: '
		print '  '+'\n  '.join(uts)
		print
	#
	if doplot:
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_pdf import PdfPages
		if decalsstyle:
			if smallplot:
				fig = plt.figure(figsize=(5,6))
				plt.subplots_adjust(0.11,0.08,0.98,0.98,0.0,0.0)
				sz1,sz2,sz3,fsz = 7,5,20,11
			else:
				fig = plt.figure(figsize=(8,10))
				plt.subplots_adjust(0.07,0.05,0.98,0.98,0.0,0.0)
				sz1,sz2,sz3,fsz = 10,12,20,12
			for _pass in range(1,4):
				ax = plt.subplot(3,1,_pass)
				grsum = tileCov[:,0,_pass-1].astype(np.int) + \
				        2*(tileCov[:,1,_pass-1].astype(np.int))
				ii = np.where(grsum==0)[0]
				plt.scatter(tiledb['TRA'][ii],tiledb['TDEC'][ii],
				            marker='+',s=sz1,c='0.7')
				ii = np.where(grsum>0)[0]
				plt.scatter(tiledb['TRA'][ii],tiledb['TDEC'][ii],
				            marker='s',
				            c=np.choose(grsum[ii],['0.5','b','y','g']),
				            edgecolor='none',s=sz2)
				if _pass==3:
					for c,lbl in zip('byg',['g','r','g+r']):
						plt.scatter(-99,-99,marker='s',s=sz3,c=c,label=lbl,
						            edgecolor='None')
					plt.legend(scatterpoints=1,ncol=3,fontsize=fsz,
					           handletextpad=0,columnspacing=1,
					           loc='upper center')
				plt.xlim(85,305)
				plt.ylim(29,62)
				if _pass==3:
					plt.xlabel('RA')
				else:
					ax.xaxis.set_ticklabels([])
				if _pass==2:
					plt.ylabel('Dec')
				plt.text(270,55,'pass %d'%_pass)
			if saveplot is not None:
				plt.savefig(saveplot)
			else:
				plt.show()
		elif byfilter:
			if brightcolors:
				passcolors = ['0.5','DarkCyan','blue','red']
			else:
				passcolors = ['0.5','#5f7ae3','#3d58c1','#031e84']
			fig = plt.figure(figsize=(12,6))
			plt.subplots_adjust(0.06,0.08,0.98,0.95,0.0,0.0)
			sz1,sz2 = 15,25
			for bi,filt in enumerate('gr'):
				ax = plt.subplot(1,2,bi+1)
				npass = tileCov[:,bi,:].sum(axis=-1)
				ii = np.where(npass==0)[0]
				plt.scatter(tiledb['TRA'][ii],tiledb['TDEC'][ii],
				            marker='+',s=sz1,c='0.7')
				ii = np.where(npass>0)[0]
				plt.scatter(tiledb['TRA'][ii],tiledb['TDEC'][ii],
				       marker='s',
				       c=np.choose(npass[ii],passcolors),
				       edgecolor='none',s=sz2)
				plt.xlim(85,305)
				plt.ylim(29,75)
				plt.title('%s band'%filt,size=14)
				if filt=='g':
					plt.ylabel('Dec',size=14)
				else:
					ax.yaxis.set_ticklabels([])
			plt.figtext(0.5,0.01,'R.A.',ha='center',size=14)
			if saveplot is not None:
				plt.savefig(saveplot)
			else:
				plt.show()
		else:
			if saveplot:
				pdf = PdfPages('bass_coverage_%s%s.pdf'%(filt,pltsfx))
			for j in range(3):
				fig = plt.figure(figsize=(10,6))
				plt.subplots_adjust(0.03,0.05,0.98,0.95)
				sz = 5 if saveplot else 20
				plt.scatter(tiledb['TRA']/15,tiledb['TDEC'],marker='s',
				            c=np.choose(nobs[:,j],
				                    ['0.9','c','DarkCyan','b','purple','m']),
				            edgecolor='none',s=sz)
				plt.plot([13+5./6,14+45./60,14+45./60,13+5./6,13+5./6],
				         [50.7,50.7,56.2,56.2,50.7],c='k')
				plt.plot([14.37,14.62,14.62,14.37,14.37],
				         [32.5,32.5,36.1,36.1,32.5],c='k')
				plt.xlim(20,5.9)
				plt.ylim(29.7,58)
				plt.title('filter %s pass %d total %d unique %d repeats %d' %
				          (filt,j+1,np.sum(nobs[:,j]),np.sum(nobs[:,j]>0),
				           np.sum(nobs[:,j]>1)))
				if saveplot:
					pdf.savefig(fig,orientation='landscape')
			if saveplot:
				pdf.close()
	return nobs,tileList

def map_to_decam_obstatus(which='good',newest=True):
	tiledb = load_tiledb()
	obsdb = load_obsdb(get_obsdb_filename(which,newest))
	#
	tid = np.array([int(tid) for tid in tiledb['TID']]).astype(np.int32)
	tabs = []
	for p in [1,2,3]:
		t = Table(dict(TILEID=tid))
		t['RA'] = tiledb['DRA%d'%p]
		t['DEC'] = tiledb['DDEC%d'%p]
		t['PASS'] = np.int16(p)
		t['EBV_MED'] = tiledb['EBV']
		tabs.append(t)
	obstatus = vstack(tabs)
	#
	obstatus['G_DONE'] = np.int16(0)
	obstatus['R_DONE'] = np.int16(0)
	obstatus['G_DATE'] = '          '
	obstatus['G_DEPTH'] = np.float32(0)
	obstatus['R_DATE'] = '          '
	obstatus['G_EXPNUM'] = np.int32(0)
	obstatus['R_EXPNUM'] = np.int32(0)
	obstatus['R_DEPTH'] = np.float32(0)
	#
	obsdb = Table(obsdb).group_by(['tileId','ditherId','filter'])
	for k,g in zip(obsdb.groups.keys,obsdb.groups):
		tileId,ditherId,filt = k
		if tileId not in tid:
			# not a BASS tile
			continue
		j = np.where((obstatus['TILEID']==tileId) & 
		             (obstatus['PASS']==ditherId))[0][0]
		fn = g['fileName'][-1]
		expNum = np.int32(fn[1:5]+fn[6:])
		utDate = Time(g['mjd'][-1],format='mjd').iso[:10]
		if filt=='g':
			obstatus['G_DONE'][j] = 1
			obstatus['G_DATE'][j] = utDate
			obstatus['G_EXPNUM'][j] = expNum
		elif filt=='r':
			obstatus['R_DONE'][j] = 1
			obstatus['R_DATE'][j] = utDate
			obstatus['R_EXPNUM'][j] = expNum
	obstatus.write('bass-tiles_obstatus.fits',overwrite=True)

def load_etc_results_file(resultsf='result.txt'):
	names = ['fileName','ra','dec','magLimNie','expTime','finalCal',
	         'skyFlux','skyRms','seeing','airmass','filter',
	         'obsCal','obsSkyFlux','obsSkyRms','obsSeeing','obsAirmass',
	         'E(B-V)','calcExpTime',
	         'magLimPoint1','magLimExt1','magLimPoint3','magLimExt3']
	return Table.read(resultsf,format='ascii',names=names)

def load_etc_maglim_file(maglimf='maglim_result_2016.txt',obsdb=None):
	if obsdb is None:
		obsdb = Table(load_obsdb())
	names = ['fileName','expTime','airmass','obsCal','seeing',
	         'skyFlux','skyRms','E(B-V)',
	         'magLimETC','magLimOBS','magLimOBSCorr']
	tab = Table.read(maglimf,format='ascii.basic',names=names)
	tab['fileName'] = [ reform_filename(fn) for fn in tab['fileName'] ]
	tab.remove_column('expTime') 
	tab = join(obsdb,tab,join_type='left',keys='fileName')
#	ii = files2tiles(obsdb,tab['fileName'])
#	print '%d tiles total, %d images missing from tile db' % \
#	            (len(ii),np.sum(ii<0))
#	m = np.where(ii>=0)[0]
#	#print tab[ii<0]
#	ii = ii[m]
#	tab = tab[m]
	# add some fields
	tab['zeroPointETC'] = 25 - 2.5*np.log10(tab['expTime']) + tab['obsCal']
#	nominalZp = np.choose(obsdb['filter'][ii]=='g',[25.38,25.55])
	nominalZp = np.choose(obsdb['filter']=='g',[25.38,25.55])
	tab['transparencyETC'] = nominalZp - tab['zeroPointETC']
	tab['skyAduPerSec'] = tab['skyFlux'] / tab['expTime']
#	tab.remove_columns(['fileName','expTime']) # these are repeats
#	return hstack([tab,Table(obsdb)[ii]])
	return tab

def cfhtw3_tiles(observed=True):
	w3west,w3east = 15*(13.+50/60.), 15*(14+45./60)
	w3south,w3north = 50.7, 56.2
	return region_tiles(w3west,w3east,w3south,w3north,observed=observed)

def ndwfs_tiles(observed=True):
	ndwest,ndeast = 15*14.37, 15*14.62
	ndsouth,ndnorth = 32.5, 36.1
	return region_tiles(ndwest,ndeast,ndsouth,ndnorth,observed=observed)

def panstarrs_md_tiles(mdfield,observed=True):
	mdfields = {'MD03':(130.592,+44.317),
                'MD05':(161.917,+58.083),
                'MD06':(185.000,+47.117),
                'MD07':(213.704,+53.083),
                'MD08':(242.787,+54.950)}
	ra,dec = mdfields[mdfield]
	dra = 3.5/np.cos(np.radians(dec))
	return region_tiles(ra-dra,ra+dra,dec-3.6,dec+3.6,
	                    observed=observed)

def deep2_tiles(deep2f,observed=True):
	d2f = {'Deep2F1':(214.25,52.5),
	       'Deep2F2':(253.0 ,34.917),
	       'Deep2F3':(352.5 ,0.0),
	       'Deep2F4':( 37.5 ,0.0)}
	ra,dec = d2f[deep2f]
	dra = 1.5/np.cos(np.radians(dec))
	return region_tiles(ra-dra,ra+dra,dec-1.6,dec+1.6,
	                    observed=observed)

def checkarchive(logf,archivef,outfn=None):
	if True:
		bassfiles = []
		with open(archivef) as basslog:
			for l in basslog:
				utdfn = re.match('.* /data/primefocus/bass/(.*)\/(.*)\n',l).groups()
				bassfiles.append(utdfn)
				bassfiles.append(utdfn)
		_archive = Table(rows=bassfiles,names=('utDir','DTACQNAM'))
		archive = unique(_archive,keys='DTACQNAM')
		print 'archive log returned %d unique files out of %d' % \
		            (len(archive),len(_archive))
	_log = Table.read(logf)
	log = unique(_log,keys='DTACQNAM')
	print 'NERSC log returned %d unique files out of %d' % (len(log),len(_log))
	# join the tables based on the unique (?!) file handle
	t = join(archive,log,keys='DTACQNAM',join_type='outer')
	if outfn:
		outf = open(outfn,'w')
	else:
		outf = sys.stdout
	if True:
		is2016 = [str(utd).startswith('2016') for utd in t['utDir_1']]
		initmask = np.array(is2016)
	# files in archive but not on NERSC
	not_on_nersc = t['frameIndex'].mask
	# files on NERSC but not archive (???)
	not_in_archive_log = t['utDir_1'].mask
	# summary info about missing files
	utds = np.unique(t['utDir_1'][initmask&not_on_nersc])
	utds = np.array(utds)
	for utd in utds:
		if 'bad' in utd: continue # files moved manually
		ii = np.where(t['utDir_1']==utd)[0]
		missing = not_on_nersc[ii]
		print '%s missing %d/%d files' % (utd,missing.sum(),len(ii))
		for i in ii[missing]:
			outf.write('%s %s\n' % (t['utDir_1'][i],t['DTACQNAM'][i]))
	if outfn:
		outf.close()

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--obsdb",action="store_true",
	                    help="build the observations database")
	parser.add_argument("--newest",action="store_true",
	                    help="include most recent observations")
	parser.add_argument("-t","--tiles",type=str,default="good",
	                    help="which set of tiles to include (good|all)")
	parser.add_argument("--update",action="store_true",
	                    help="update the tile list by remote download")
	parser.add_argument("--summary",action="store_true",
	                    help="print summary information")
	parser.add_argument("-p","--plot",action="store_true",
	                    help="also produce summary plot")
	parser.add_argument("--smallplot",action="store_true",
	                    help="small summary plot")
	parser.add_argument("-m","--mjd",type=str,
	                    help="MJD range (mjd1,mjd2), '*' for any")
	parser.add_argument("--start2016",action="store_true",
	                    help="start with Nov 1 2015 ('2016' data)")
	parser.add_argument("-u","--utdate",type=str,
	                    help="UT range (ut1,ut2), '*' for any")
	parser.add_argument("-d","--date",type=str,
	                    help="local date range (night1,night2), '*' for any")
	parser.add_argument("--plotfile",type=str,
	                    help="filename to save plot in")
	parser.add_argument("--legacyfield",type=str,
	                    help="select tiles within legacy field")
	parser.add_argument("--byfilter",action="store_true",
	                    help="plot coverage by filter instead of by pass")
	parser.add_argument("--bright",action="store_true",
	                    help="plot using higher contrast colors")
	parser.add_argument("-v","--verbose",action="store_true",
	                    help="increase verbosity")
	args = parser.parse_args()
	mjds = [None,None]
	if args.mjd is not None:
		mjds = [ int(d) 
		            if d!='*' else None for d in args.mjd.split(',') ]
	elif args.utdate is not None:
		mjds = [ int(Time(d).mjd) 
		            if d!='*' else None for d in args.utdate.split(',') ]
	elif args.date is not None:
		mjds = [ int((Time(d)+TimeDelta(1,format='jd')).mjd) 
		            if d!='*' else None for d in args.date.split(',') ]
	if args.start2016:
		mjds[0] = 57327
	if len(mjds)==1:
		mjds = mjds*2
	if args.obsdb:
		build_obsdb(update=args.update,which=args.tiles,newest=args.newest)
	elif args.summary:
		if args.legacyfield:
			if args.legacyfield.upper() == 'CFHTW3':
				tiles = cfhtw3_tiles(observed=False)
			elif args.legacyfield.upper() == 'NDWFS':
				tiles = ndwfs_tiles(observed=False)
			elif args.legacyfield.upper().startswith('MD'):
				tiles = panstarrs_md_tiles(args.legacyfield.upper(),
				                           observed=False)
			elif args.legacyfield.startswith('Deep2'):
				tiles = deep2_tiles(args.legacyfield,observed=False)
			else:
				raise ValueError
			tiles = np.array([int(tid) for tid in tiles['TID']])
			print '\nnumber of tiles covering %s: %d' % \
			           (args.legacyfield,len(tiles))
		else:
			tiles = None
		obs_summary(which=args.tiles,newest=args.newest,tiles=tiles,
		            mjdstart=mjds[0],mjdend=mjds[1],
		            doplot=args.plot,smallplot=args.smallplot,
		            saveplot=args.plotfile,decalsstyle=not args.byfilter,
		            byfilter=args.byfilter,brightcolors=args.bright,
		            verbose=args.verbose)
	#kwargs = {} if len(sys.argv)==1 else {'dirs':sys.argv[1]}
	#print kwargs
	#nersc_archive_list(**kwargs)


