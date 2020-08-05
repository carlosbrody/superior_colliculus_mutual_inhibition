#!/usr/bin/env bash

#SBATCH -J 'spockFarm'
#SBATCH -o  log-spockFarm-%j.out
#SBATCH -p Brody
#SBATCH --time 96:00:00
#SBATCH -c 1

# load the julia module
module load julia/1.2.0

# print out some info for the log file
echo "Slurm Job ID, unique: $SLURM_JOB_ID"
echo "Slurm Array Task ID, relative: $SLURM_ARRAY_TASK_ID"

# call the julia script which will start the farm
julia $1 --bspockID ${SLURM_ARRAY_JOB_ID}_$SLURM_ARRAY_TASK_ID --sleepFor 10 --onSpock --hostname spock

# move the log file to an input-dependent location
# mv log-array-example-${SLURM_JOB_ID}.out $1/log-array-example-${SLURM_JOB_ID}.out
