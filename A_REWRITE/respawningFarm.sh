#! /bin/csh

# Runs a julia file n times in the background.  Useful for starting a bunch of identical julias going
# For example:
#
# unix>  ./respawningFarm.sh   julia_file.jl   60
#
# Will start 60 identical processes in the background.


if ( $#argv != 2 ) then
    echo Sorry, need two arguments, the name of the julia file to farm, and the number of farms to start
    exit
endif

setenv HOST `hostname`

set j=1
while ( $j <= $argv[2] )
   echo "Will run: rm julia_outs_${HOST}_$j ; julia $argv[1] $j $argv[2] >> julia_outs_${HOST}_$j &"
   rm julia_outs_${HOST}_$j
   julia $argv[1] $j $argv[2] >> julia_outs_${HOST}_$j &
   sleep 3
   @ j++
end

while ( 2 > 1 )  # permanently
   set j=1
   while ( $j <= $argv[2] )  # go through all started processes
      set myproc=`ps aux | grep "julia $argv[1] $j $argv[2]" | grep -v grep | wc -l`
      if ( $myproc == 0 ) then  # if we couldn't find one
         echo "    Respawning :  julia $argv[1] $j $argv[2] >> julia_outs_${HOST}_$j &"
         julia $argv[1] $j $argv[2] >> julia_outs_${HOST}_$j &
         sleep 3
      endif
      @ j++
   end
   sleep 5  # pause before doing another check of all processes
end
