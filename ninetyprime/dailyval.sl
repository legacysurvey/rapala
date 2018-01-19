#!/bin/bash -l

##SBATCH -C haswell
#SBATCH -p regular
#SBATCH -t 01:00:00

##SBATCH -p debug
##SBATCH -t 00:30:00

#SBATCH -N 1 

# Generate the nightly logs for incoming raw data. The -u switch is used
# to limit which directories are scanned, this speeds things up by 
# allowing the logs to be generated for only the most recent data
cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior
srun -n 1 python detchar.py -n -d ~/bok/BOK_Raw/ -l basslogs -u "201711*" 


# the following commands were used to generate QA files only for the
# purpose of checking overheads during a full season

#srun -n 1 -c 24 python detchar.py -n -d ~/bok/BOK_Raw/ --logfile ../survey/logs2016.fits --nogainrn --nobitcheck -o bassqa2016.fits --nproc 24

#srun -n 1 -c 16 python detchar.py -n -d ~/bok/BOK_Raw/ --logfile ../survey/logs2015.fits --nogainrn --nobitcheck -o bassqa2015.fits --nproc 24 

#srun -n 1 -c 16 python detchar.py -n -d ~/bok/BOK_Raw/ --logfile ../survey/logs2015.fits --nogainrn --nobitcheck -o bassqa2015_1.fits --nproc 24 --numsplit 2 --splitnum 0
#srun -n 1 -c 16 python detchar.py -n -d ~/bok/BOK_Raw/ --logfile ../survey/logs2015.fits --nogainrn --nobitcheck -o bassqa2015_2.fits --nproc 24 --numsplit 2 --splitnum 1

