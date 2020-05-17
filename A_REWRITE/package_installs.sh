#!/bin/sh

# This shell script will take a brand new Deebian machine (like the ones in Google Compute Engine)
# and will install all the packages needed to get Julia up and running for parameter minimization.
#

# ON UBUNTU 16.04:
sudo apt-get -y install build-essential
sudo apt-get update
sudo apt-get -y upgrade
sudo apt-get -y install gfortran
sudo apt-get -y install cmake
sudo apt-get -y install m4
sudo apt-get -y install git-core
#
#  PUT YOUR OWN NAME AND EMAIL HERE !!!!
#
git config --global user.email "carlosbrody@gmail.com"
git config --global user.name "Carlos Brody"
sudo apt-get -y install emacs
sudo apt-get -y install python-dev
sudo apt-get -y install python-matplotlib
sudo apt-get -y install csh
sudo apt-get -y install hdf5-tools

# wget https://julialang-s3.julialang.org/bin/linux/x64/0.5/julia-0.5.2-linux-x86_64.tar.gz
# tar xvfz julia-0.5.2-linux-x86_64.tar.gz
# sudo ln -s /home/carlosbrody/julia-f4c6c9d4bb/bin/julia /usr/bin/julia.0.5.2

# wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.1-linux-x86_64.tar.gz
# tar xvfz julia-0.6.1-linux-x86_64.tar.gz
# sudo ln -s /home/carlosbrody/julia-0d7248e2ff/bin/julia /usr/bin/julia


#
# CUDA Stuff, only if a GPU is on
#

if [ -f  /usr/bin/lspci ] ;  then
    if ( lspci | grep -qi nvidia ) ; then
	wget http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/cuda-repo-ubuntu1604_9.1.85-1_amd64.deb
	sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
	sudo dpkg -i cuda-repo-ubuntu1604_9.1.85-1_amd64.deb
	sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
	sudo apt-get update
	sudo apt-get -y install cuda

	# GET AND BUILD JULIA
	git clone git://github.com/JuliaLang/julia.git
	cd julia
	git checkout v0.6.2
	make
	cd ..
	sudo ln -s /home/carlosbrody/julia/julia /usr/bin/julia
    fi
else
    wget wget https://julialang-s3.julialang.org/bin/linux/x64/1.0/julia-1.0.5-linux-x86_64.tar.gz
    tar xvfz julia-1.0.5-linux-x86_64.tar.gz
    sudo ln -s /home/carlosbrody/julia-1.0.5/bin/julia /usr/bin/julia
fi




#
#  PUT YOUR OWN GITHUB USERNAME AND PASSWORD HERE !!!!
#
git clone "http://carlosbrody:gungaDID90&@github.com/carlosbrody/superior_colliculus_mutual_inhibition.git"


# Final step is to run julia to install its packages:
julia julia_start.jl

# Grow the ForwardDiff chunks *after* ForwardDiff is installed:
./grow_ForwardDiffChunks.sh
