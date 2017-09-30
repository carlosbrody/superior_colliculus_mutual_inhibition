#!/bin/csh

cat ~/.julia/v0.5/ForwardDiff/src/ForwardDiff.jl | sed 's/const CHUNK_THRESHOLD = 10/const CHUNK_THRESHOLD = 15/' > temp
mv temp ~/.julia/v0.5/ForwardDiff/src/ForwardDiff.jl


