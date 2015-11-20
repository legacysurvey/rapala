#!/usr/bin/env/python

import sys

from bokpipe.bokastrom import scamp_solve

scamp_solve(*sys.argv[1:],overwrite=True)

