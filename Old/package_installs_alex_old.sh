#!/bin/sh

# This shell script will take a brand new Deebian machine (like the ones in Google Compute Engine)
# and will install all the packages needed to get Julia up and running for parameter minimization.
#

# ON DEBIAN:
sudo apt-get -y install build-essential
sudo apt-get update
sudo apt-get -y upgrade
sudo apt-get -y install git-core
git config --global user.email "alexpiet@gmail.com"
git config --global user.name "alexpiet"
sudo apt-get -y install emacs
sudo apt-get -y install python-dev
sudo apt-get -y install csh

wget https://julialang-s3.julialang.org/bin/linux/x64/0.5/julia-0.5.2-linux-x86_64.tar.gz
tar xvfz julia-0.5.2-linux-x86_64.tar.gz 
sudo ln -sf /home/alex.piet/julia-f4c6c9d4bb/bin/julia /usr/bin/julia

sudo apt-get -y install hdf5-tools

git clone "http://alexpiet@github.com/carlosbrody/superior_colliculus_mutual_inhibition.git"

# THEN START JULIA AND RUN THIS FILE:

julia julia_start.jl




