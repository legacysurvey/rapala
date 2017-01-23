#!/usr/bin/make

.PHONY: nov15initproc nov15badpix nov15flats nov15finalproc

TOOLS := $(HOME)/dev/rapala/bokpipe/tools

#NOV15LOG := config/nov2015.fits
# used topcat to flag images to ignore (added 'good' column)
NOV15LOG := config/nov2015_mod.fits

HOST=$(shell hostname)
ifeq ($(HOST),clarno)
	NOV15RAW := $(BASSDATA)/Nov2015
else
	NOV15RAW := $(BASSDATA)/BOK_Raw
endif

NPROC := -p 6

nov15all: nov15initproc nov15badpix nov15flats nov15finalproc

nov15db:
	python $(TOOLS)/bokmkobsdb.py -o config/nov2015.fits \
	       $(BASSDATA)/Nov2015/2015111?

nov15initproc:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s oscan,bias2d,flat2d --nousepixflat \
	                    -vvv $(NPROC)

nov15badpix:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s bpmask -u 20151115 -b g -vvv

# exclude the linearity sequences at fixed positions
NOV14SKYFLATS := "33-38,45-50,57-110"

nov15flats:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s proc1 --prockey TMPPROC \
	                    --norampcorr --nobiascorr --nosavegain \
	                    --noweightmap --tmpdirout $(NPROC)
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --makeillumcorr \
	                    --tmpdirin -vvv $(NPROC)
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    -s proc2 \
	                    --nodarkskycorr --noweightmap \
	                    --skymethod polynomial --skyorder 1 \
	                    --tmpdirin --tmpdirout -vvv $(NPROC)
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14SKYFLATS) \
	                    -s skyflat \
	                    --tmpdirin --tmpdirout -vvv  $(NPROC)

# Just the science frames in validation fields
NOV14FRAMES := "33-69,98-110"

nov15finalproc:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s proc1 --norampcorr --nobiascorr $(NPROC)
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s proc2 $(NPROC)
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s wcs,cat -vv

nov15finish:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) --finish

nov15proc:
	python basschute.py --obsdb $(NOV15LOG) -r $(NOV15RAW) \
	                    --frames $(NOV14FRAMES) \
	                    -s oscan,proc1,proc2,wcs,cat \
	                    --norampcorr --nobiascorr \
	                    -vvv $(NPROC) $(XARGS)

nov15redo:
	make XARGS=-R nov15proc 

