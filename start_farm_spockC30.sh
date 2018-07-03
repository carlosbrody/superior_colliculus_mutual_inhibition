#!/usr/bin/env bash

#SBATCH -J 'farm-C30'
#SBATCH -o log-farmC30-%j.out
#SBATCH -p Brody
#SBATCH --time 10:00:00
#SBATCH -c 1


# Start this script with
# sbatch --array=0-9 ./start_farm_spockC30.sh 

# load the julia module
module load julia/0.6.3

# print out some info for the log file
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"

# Call the julia script which will start the farm
julia spock_opto_reduced_farmC30.jl $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID

# Move the log file
mv log-farmC30-${SLURM_JOB_ID}.out ../Reports/log-farmC30-${SLURM_JOB_ID}.out


