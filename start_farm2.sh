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

set j=1
while ( $j <= $argv[2] )
   echo "Will run:  julia $argv[1] > julia_outs_$j &"
   julia $argv[1] >& julia_outs_$j &
   sleep 2
   @ j++
end



