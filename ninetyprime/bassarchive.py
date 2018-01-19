#!/usr/bin/env python

import os
import sys
import glob
import subprocess
import md5

archivelogfn = "/home/primefocus/bass/archive/bassarchive.log"

dtskp = "bokpf@dtskp.kpno.noao.edu"

def check_file_size(year=''):
	'''List the files on the dtskp machine and check the output size in
	bytes. Return a list of files that are the nominal 90Prime size.'''
	rv = subprocess.check_output(["ssh",dtskp,
	                              "ls","-lR","bass/{}*".format(year)])
	dtsfiles = {}
	for l in rv.split('\n'):
		d = l.strip().split()
		try:
			if int(d[4]) == 134841600:
				dtsfiles[d[-1]] = 1
		except IndexError:
			pass
	return sorted(dtsfiles.keys())

def read_archive_log():
	'''Read in the log of files transferred to dtskp. Returns a dict where
	the keys are the filenames and the values are the MD5 hash strings.'''
	archivelog = {}
	with open(archivelogfn) as archivelogf:
		for l in archivelogf:
			try:
				md5h,fn = l.strip().split()
			except:
				print "BASSARCHIVE: ERROR reading ",l
				continue
			fn = os.path.basename(fn)
			if fn in archivelog and archivelog[fn] != md5h:
				if False:
					print fn,archivelog[fn],md5h
			else:
				archivelog[fn] = md5h
	return archivelog

def get_files_to_transfer(currentdirs,archivelog,dtsfiles):
	'''List files in the Bok raw data directory and check them against
	the archive log file. Return a list of files in the raw data directory
	but not in the log; i.e., the files that haven't been transferred yet.'''
	xferfiles = []
	for d in currentdirs:
		currentfiles = sorted(glob.glob(os.path.join(d,"d????.????.fits")))
		for f in currentfiles:
			fn = os.path.basename(f)
			if not fn in archivelog:
				xferfiles.append(f)
			if dtsfiles is not None and fn not in dtsfiles:
				print "BASSARCHIVE: {} not on dtskp".format(fn)
	return xferfiles

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-d","--datadir",type=str,
	                    default="/data/primefocus/bass/",
	                    help="top-level data directory")
	parser.add_argument("-y","--year",type=str,
	                    default="2018",
	                    help="restrict archive check to single year")
	parser.add_argument("--lastnight",action='store_true',
	                    help="only run archiving for last night of data")
	parser.add_argument("--checksize",action='store_true',
	                    help="additionally check file size on dtskp")
	parser.add_argument("--nocopy",action='store_true',
	                    help="don't copy any files, just print output")
	parser.add_argument("--scp",action='store_true',
	                    help="use scp instead of rsync")
	args = parser.parse_args()

	# read in the list of files in the archive log
	archivelog = read_archive_log()

	# get the list of nightly directories in the raw data directory
	dirstr = os.path.join(args.datadir,"{}*".format(args.year))
	currentdirs = sorted(glob.glob(dirstr))

	# if requested, check the file sizes of files on dtskp and return a
	# list of those that have the correct size (catches interruped transfers)
	if args.checksize:
		dtsfiles = check_file_size(args.year)
	else:
		dtsfiles = None 

	# get the list of files that need to be copied to dtskp
	xferfiles = get_files_to_transfer(currentdirs,archivelog,dtsfiles)

	# if requested, only do the last night's data
	if args.lastnight:
		# only the last night
		currentdirs = [currentdirs[-1]]

	try:
		dtsutdirs = subprocess.check_output(["ssh",dtskp,
		                                     "ls","-1d",
		                                     "bass/{}*".format(year)])
	except:
		dtsutdirs = ''
	dtsutdirs = [os.path.basename(d) for d in dtsutdirs.split("\n") ]

	if len(xferfiles) == 0:
		print 'BASSARCHIVE: nothing to do!'
		sys.exit(0)

	# Finally do the transfer. Open the archive log files, iterate through
	# the list of files to transfer, and copy files one-by-one.
	FNULL = open(os.devnull, 'w')
	with open(archivelogfn,"a") as archivelogf:
		for f in xferfiles:
			d,fn = os.path.split(f)
			d,utdir = os.path.split(d)
			# if the night directory doesn't exist on dtskp, create it
			if utdir not in dtsutdirs:
				cmd = ["ssh",dtskp,"mkdir","bass/"+utdir]
				print "BASSARCHIVE: "+(' '.join(cmd))
				_ = subprocess.check_output(cmd)
				dtsutdirs.append(utdir)
			if args.scp:
				cmd = ["scp","-r",f,dtskp+":bass/"+utdir+"/"]
			else:
				cmd = ["rsync","-av",f,dtskp+":bass/"+utdir+"/"]
			print "BASSARCHIVE: "+(' '.join(cmd))
			if not args.nocopy:
				# copy the file and add the MD5 hash to the log file
				rv = subprocess.call(cmd,stdout=FNULL)
				rv = subprocess.check_output(["md5sum",f]).strip()
				archivelogf.write(rv+"\n")

