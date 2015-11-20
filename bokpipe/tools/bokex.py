#!/usr/bin/env/python

import sys

from bokpipe.bokextract import sextract

sextract(*sys.argv[1:],overwrite=True)

