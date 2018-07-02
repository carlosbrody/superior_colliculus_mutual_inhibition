#!/usr/bin/env bash

#SBATCH -J 'setup_julia'
#SBATCH -o  log_setup_julia-%j.out
#SBATCH -p Brody
#SBATCH --time 00:30:00
#SBATCH --mem 1000
#SBATCH -c 1

module load julia/0.6.3

julia setup.jl

