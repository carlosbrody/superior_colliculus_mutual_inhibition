#!/usr/bin/env bash

#SBATCH -J 'julia-example'
#SBATCH -o  log-julia-example-%j.out
#SBATCH -p Brody
#SBATCH --time 00:10:00
#SBATCH --mem 1000
#SBATCH -c 1

module load julia/0.6.3
module load

julia example.jl

