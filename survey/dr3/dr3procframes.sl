#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1 
#SBATCH -t 02:00:00

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

#srun -n 1 -c 64 python bass_ccds.py -p --dr3 -m 64 -x wcs_failed bass_Nov2015toFeb2016.fits

# pick up the WCS failures with a larger (and slower) search radius
srun -n 1 -c 24 python bass_ccds.py -p --dr3 -m 24 bass_Nov2015toFeb2016.fits --searchradius 1.0 

#srun -n 1 -c 64 python bass_ccds.py -p --dr3 -m 64 -x wcs_failed bass_Nov2015toFeb2016.fits -n 1,1000

