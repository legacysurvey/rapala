#!/usr/bin/make

TOOLS := $(HOME)/dev/rapala/bokpipe/tools

#NOV15LOG := config/nov2015.fits
# used topcat to flag images to ignore (added 'good' column)
NOV15LOG := config/nov2015_mod.fits
NOV15RAW := $(BASSDATA)/Nov2015

nov15db:
	python $(TOOLS)/bokmkobsdb.py -o config/nov2015.fits \
	       $(BASSDATA)/Nov2015/2015111?

nov15proc:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s oscan,bias2d,flat2d -vvv #-p 4 

nov15badpix:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s bpmask -u 20151115 -b g -vvv

nov15flats:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s proc1 --prockey TMPPROC \
	                    --norampcorr --nosavegain --tmpdirout #-p 4
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --tmpdirin --makeillumcorr -vvv
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s proc2 \
	                    --nodarkskycorr --skymethod polynomial --skyorder 1 \
	                    --tmpdirin --tmpdirout -vvv #-p 4
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s skyflat \
	                    --tmpdirin --tmpdirout -vvv #-p 4

# Just the science frames in validation fields
NOV14FRAMES := "33-69,98-110"

nov15final:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s proc1 --norampcorr
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s proc2 -b g --nodarkskycorr
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s proc2 -b bokr --nodarkskycorr
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s wcs,cat -vv

nov15finish:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) --finish

