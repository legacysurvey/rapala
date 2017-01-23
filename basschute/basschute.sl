#!/bin/bash -l

#SBATCH -p regular
#SBATCH -t 05:00:00
#SBATCH -J basschute
#SBATCH -N 1 

##SBATCH -p shared
##SBATCH -t 03:00:00
##SBATCH -J basschute
###SBATCH -n 1 
##SBATCH --mem=10GB

##SBATCH -p debug
##SBATCH -t 00:30:00
##SBATCH -N 1 

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior
##srun -n 1 make

##srun -n 1 make nov15proc
##srun -n 1 make nov15badpix
##srun -n 1 make nov15flats

export UTDATE=20160703 

if [ "$NERSC_HOST" == "edison" ]
then
	export NPROC=24
	export MAXMEM=50
	export BASSRDXDIR=.
fi

if [ "$NERSC_HOST" == "cori" ]
then
	export NPROC=32
	#export NPROC=10
	export MAXMEM=100
fi

echo "redux dir is $BASSRDXDIR"

export XARGS="-p $NPROC --maxmem $MAXMEM"
# include cals
#export XARGS="-p $NPROC --maxmem $MAXMEM --frames 16601-16640"
# only objects
#export XARGS="-p $NPROC --maxmem $MAXMEM --frames 16621-16640"

srun -n 1 -c $NPROC make -f Makefile.ut all

#srun -n 1 -c $NPROC make -f Makefile.ut initproc badpix
#srun -n 1 -c $NPROC make -f Makefile.ut flats

