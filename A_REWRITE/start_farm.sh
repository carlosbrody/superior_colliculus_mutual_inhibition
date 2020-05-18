#! /bin/csh

# Runs a julia file n times in the background.  Useful for starting a bunch of identical julias going
# For example:
#
# unix>  ./start_farm.sh   julia_file.jl   60
#
# Will start 60 identical processes in the background.

if ( $#argv != 2 ) then
    echo Sorry, need two arguments, the name of the julia file to farm, and the number of farms to start
    exit
endif

setenv HOST `hostname`

set j=1
while ( $j <= $argv[2] )
   echo "Will run:  julia $argv[1] $j $argv[2] > julia_outs_${HOST}_$j &"
   # julia $argv[1] $j $argv[2] > julia_outs_$HOST_$j &
   sleep 3
   @ j++
end
