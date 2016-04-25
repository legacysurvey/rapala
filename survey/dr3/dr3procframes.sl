#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1 
#SBATCH -t 02:00:00

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

srun -n 1 -c 64 python bass_ccds.py -p --dr3 -m 16 -x wcs_failed bass_Nov2015toFeb2016.fits

#srun -n 1 -c 64 python bass_ccds.py -p --dr3 -m 16 -x wcs_failed bass_Nov2015toFeb2016.fits -n 1,1000

