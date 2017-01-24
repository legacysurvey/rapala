#!/usr/bin/env python

from astropy.table import Table

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("obsdb",type=str,
                    help="input observations database")
args = parser.parse_args()

fields = ['frameIndex','fileName','imType',
          'filter','objName','expTime']

obsdb = Table.read(args.obsdb)
obsdb[fields].more()

