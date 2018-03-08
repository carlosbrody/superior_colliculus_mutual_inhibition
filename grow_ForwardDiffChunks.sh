#!/bin/csh
#
# Run on its own, this script will check for a ForwardDiff on Julia 0.5 and on a Julia 0.6 installation
# in the user's home directory, and will 
# update them to have default chunk sizes of 15 as opposed to 10.
#

set J5_file = ~/.julia/v0.5/ForwardDiff/src/ForwardDiff.jl
set J6_file = ~/.julia/v0.6/ForwardDiff/src/prelude.jl

if ( -e $J5_file ) then
    set tempfile = /tmp/grow_chunk_temp.jl
    cat $J5_file | sed 's/const CHUNK_THRESHOLD = 10/const CHUNK_THRESHOLD = 20/' > $tempfile
    mv $tempfile $J5_file
endif


if ( -e $J6_file ) then
    set tempfile = /tmp/grow_chunk_temp.jl
    cat $J6_file | sed 's/const DEFAULT_CHUNK_THRESHOLD = 10/const DEFAULT_CHUNK_THRESHOLD = 20/' > $tempfile
    mv $tempfile $J6_file
endif



