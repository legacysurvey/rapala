#!/usr/bin/env python

import sys
import bokutil

infile,op,operand,outfile = sys.argv[1:]

imarith = bokutil.BokImArith(op,operand,output_map=lambda f: outfile)

imarith.process_files([infile])


