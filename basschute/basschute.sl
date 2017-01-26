#!/bin/bash -l

#SBATCH -p regular
#SBATCH -t 01:30:00
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

if [ "$NERSC_HOST" == "edison" ]
then
	if [[ -z "$NPROC" ]]; then
		export NPROC=24
	fi
	export MAXMEM=50
fi

if [ "$NERSC_HOST" == "cori" ]
then
	if [[ -z "$NPROC" ]]; then
		export NPROC=32
	fi
	export MAXMEM=100
fi

echo "redux dir is $BASSRDXDIR"

export XARGS="--maxmem $MAXMEM $XARGS"

srun -n 1 -c $NPROC make quickproc

#srun -n 1 -c $NPROC make all_detrend

#srun -n 1 -c $NPROC make initproc
#srun -n 1 -c $NPROC make badpix
#srun -n 1 -c $NPROC make proc1 
#srun -n 1 -c $NPROC make makeillum
##srun -n 1 -c $NPROC make flats
#srun -n 1 -c $NPROC make proc2

