#!/usr/bin/make

.PHONY: initproc badpix flats proc redo

ifndef DATE
	$(error DATE must be set!)
endif
UTARGS := --night $(DATE)

ifdef LOGFILE
	LOGARGS := --obsdb $(LOGFILE)
endif

ifndef RAWDATA
	RAWDATA := $(BASSDATA)/BOK_Raw
endif
DATAARGS := -r $(RAWDATA)

ifdef NPROC
	MPARGS = -p $(NPROC)
endif

ifndef VERBOSE
	VERBOSE := -v
endif

ifdef BANDS
	BANDARGS := -b $(BANDS)
endif

ifdef BASS2015
	PROCARGS := --nousepixflat --rampcorr --noweightmap
else
	PROCARGS := --nousepixflat --fixsaturation --nobiascorr --noweightmap
endif

ifndef DOSKYFRINGE
	PROCARGS := $(PROCARGS) --noskyflatcorr --nofringecorr
endif

INITARGS := $(LOGARGS) $(DATAARGS) $(UTARGS) $(BANDARGS) \
            $(MPARGS) $(VERBOSE) 

all_detrend: initproc badpix proc1 makeillum flats proc2

obsdb:
	python $(BOKTOOLS)/bokmkobsdb.py -o $(LOGFILE) $(RAWDATA)/$(DATE)

# Overscan-subtract all images, generate 2D biases and dome flats
initproc:
	python basschute.py $(INITARGS) $(PROCARGS) --badpixdir $(BASSDIR) \
	                    -s oscan,bias2d,flat2d \
	                    $(XARGS)

# Generate the bias ramp image (<=2015 data)
biasramp:
	python basschute.py $(INITARGS) \
	                    -s ramp \
	                    $(XARGS)

# XXX copy in a master bp mask from config dir
badpix:
	python basschute.py $(INITARGS) \
	                    --makebpmask DomeFlat201611060110g \
	                    $(XARGS)

# First-round processing: bias/domeflat correction, combine into CCD extensions
proc1:
	python basschute.py $(INITARGS) $(PROCARGS) \
	                    -s proc1 \
	                    $(XARGS)

# Make the illumination correction image
makeillum:
	python basschute.py $(INITARGS) \
	                    -s illum 
	                    $(XARGS)

#
# Sky flat generation (processing output to temp directory)
#

#  ... apply the illumination correction to the sky flat images
flats_illumcorr:
	python basschute.py $(INITARGS) \
	                    -s proc2 --prockey TMPPRO2 \
	                    --nofringecorr --noskyflatcorr \
	                    --noskysub --noweightmap \
	                    --tmpdirout 
	                    $(XARGS)

#  ... make fringe masters from sky flat images
flats_makefringe:
	python basschute.py $(INITARGS) \
	                    -s fringe \
	                    --tmpdirin --tmpdirout 
	                    $(XARGS)

# ... apply fringe correction to sky flat images
flats_fringeskycorr:
	python basschute.py $(INITARGS) \
	                    -s proc2 --prockey TMPPRO3 \
	                    --noillumcorr --noskyflatcorr --noweightmap \
	                    --skymethod polynomial --skyorder 1 \
	                    --tmpdirin --tmpdirout 
	                    $(XARGS)

# ... combine temp processed images to make sky flat
flats_makeskyflat:
	python basschute.py $(INITARGS) \
	                    -s skyflat \
	                    --tmpdirin --tmpdirout 
	                    $(XARGS)

# all the steps to generate sky flats in one target
flats: flats_illumcorr \
       flats_makefringe flats_fringeskycorr \
       flats_makeskyflat

# Second-round processing: apply illumination, skyflat, and fringe corrections
#  and do sky subtraction
proc2:
	python basschute.py $(INITARGS) $(PROCARGS) \
	                    -s proc2 \
	                    $(XARGS)

# Perform individual processing steps as listed in STEPS
steps:
	python basschute.py $(INITARGS) $(PROCARGS) -s $(STEPS) $(XARGS)

# Assuming cals already exist, perform all the processing steps on science ims
procall:
	python basschute.py $(INITARGS) $(PROCARGS) \
	                    -s oscan,proc1,proc2 -t object \
	                    $(XARGS)

# quick-processing that doesn't include sky flat / fringe corrections
quickproc:
	python basschute.py $(INITARGS) $(PROCARGS) --badpixdir $(BASSDIR) \
	                    -s oscan,bias2d,flat2d,proc1,illum,proc2 \
	                    $(XARGS)


# Obtain astrometric solutions
wcs:
	python basschute.py $(INITARGS) -s wcs --gaia $(XARGS)


# Generate object catalogs and PSF models with sextractor+psfex
catalogs:
	python basschute.py $(INITARGS) -s cat $(XARGS)


#
# Diagnostic stuff
#

# Make PNG images for inspection purposes
images:
	python basschute.py $(INITARGS) --images $(XARGS)

