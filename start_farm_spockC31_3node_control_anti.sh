#!/usr/bin/env bash

#SBATCH -J 'fC31a'
#SBATCH -o log-farmC31-control-anti-%j.out
#SBATCH -p Brody
#SBATCH --time 96:00:00
#SBATCH -c 1


# Start this script with
# sbatch --array=0-9 ./start_farm_spockC31.sh 

# load the julia module
module load julia/0.6.3

# print out some info for the log file
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"

# Call the julia script which will start the farm
julia farmC31_3node_control_anti.jl $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID

# Move the log file
mv log-farmC31-control-anti-${SLURM_JOB_ID}.out ../Reports/log-farmC31-control-anti-${SLURM_JOB_ID}.out


