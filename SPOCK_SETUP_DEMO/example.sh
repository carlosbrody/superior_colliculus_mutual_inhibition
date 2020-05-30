#!/usr/bin/env bash

#SBATCH -J 'julia-example'
#SBATCH -o  log-julia-example-%j.out
#SBATCH -p Brody
#SBATCH --time 00:10:00
#SBATCH --mem 1000
#SBATCH -c 1

echo "I'm alive!" >> imalive.txt

module load julia/1.2.0
module load

julia example.jl

