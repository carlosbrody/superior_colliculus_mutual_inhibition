#!/bin/bash

for i in `seq 1 64`
do
    cat julia_outs_$i | grep Came | cut -d "d" -f 2 > crashed_run_$i.jl
    if ! [ -s crashed_run_$i.jl ]
    then
        rm crashed_run_$i.jl
    fi
done
