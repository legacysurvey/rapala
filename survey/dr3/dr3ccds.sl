#!/bin/bash -l

#SBATCH -N 1 
#SBATCH -p regular
#SBATCH -t 01:30:00

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

BOK_TOOLS_DIR=$HOME/github/rapala/bokpipe/tools
BASS_DATA_DIR=$BASSDATA/BOK_Raw
#BASS_DR3_DIR=$BASSDATA/ccds_files
BASS_DR3_DIR=./
BASS_DR3_FRAMES=$BASS_DR3_DIR/bass_Nov2015toFeb2016.fits

srun -n 1 python $BOK_TOOLS_DIR/bokmkobsdb.py $BASS_DATA_DIR/201511* $BASS_DATA_DIR/20160[12]* -o $BASS_DR3_FRAMES -e "DTACQNAM:S15"

# use quick-processed raw frames
#srun -n 1 python bass_ccds.py $BASS_DR3_FRAMES -r -f $BASS_DR3_DIR/bass-ccds-dr3raw.fits

# use NAOC-reduced images
srun -n 1 python bass_ccds.py $BASS_DR3_FRAMES -f $BASS_DR3_DIR/bass-ccds-dr3proc.fits

