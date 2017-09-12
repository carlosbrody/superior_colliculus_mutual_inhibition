#! /bin/csh

if ( $#argv != 2 ) then
    echo Sorry, need two arguments, the name of the julia file to farm, and the number of farms to start
    exit
endif

set j=1
while ( $j <= $argv[2] )
   echo "Will run:  julia $argv[1] > julia_outs_$j &"
   julia $argv[1] > julia_outs_$j &
   sleep 3
   @ j++
end



