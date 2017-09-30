#!/bin/sh

# This shell script will take a brand new Deebian machine (like the ones in Google Compute Engine)
# and will install all the packages needed to get Julia up and running for parameter minimization.
#

# ON DEBIAN:
sudo apt-get -y install build-essential
sudo apt-get update
sudo apt-get -y upgrade
sudo apt-get -y install git-core
#
#  PUT YOUR OWN NAME AND EMAIL HERE !!!!
#
git config --global user.email "YOUR@EMAIL"
git config --global user.name "YOUR NAME"
sudo apt-get -y install emacs
sudo apt-get -y install python-dev
sudo apt-get -y install python-matplotlib
sudo apt-get -y install csh
sudo apt-get -y install hdf5-tools

wget https://julialang-s3.julialang.org/bin/linux/x64/0.5/julia-0.5.2-linux-x86_64.tar.gz
tar xvfz julia-0.5.2-linux-x86_64.tar.gz 
sudo ln -s /home/carlosbrody/julia-f4c6c9d4bb/bin/julia /usr/bin/julia


#
#  PUT YOUR OWN GITHUB USERNAME AND PASSWORD HERE !!!!
#
git clone "http://your_github_usernam:your_github_password@github.com/carlosbrody/superior_colliculus_mutual_inhibition.git"


./grow_ForwardDiffChunks.sh

# Final step is to run julia to install its packages:
julia julia_start.jl




