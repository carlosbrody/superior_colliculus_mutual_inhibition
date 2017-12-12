#! /bin/csh

# This script goes kills -9 ALL processes matching a name. Crazy dangerous!

if ( $#argv != 1 ) then
    echo Sorry, need one arguments, the name of the process(es) to be killed
    echo For example:  ./kill_by_name.sh julia
    exit
endif

ps -ef | grep $argv[1] | grep -v grep | grep -v kill_by_name | awk '{print $2}' | xargs kill -9
