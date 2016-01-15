#!/usr/bin/env python

import sys
from bokpipe import bokproc

infile,op,operand,outfile = sys.argv[1:]

imarith = bokproc.BokImArith(op,operand,output_map=lambda f: outfile)

imarith.process_files([infile])


