# Run this file in julia, and it will add and build all the necessary packages for cost function minimization

ENV["PYTHON"]="/usr/bin/python"
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
Pkg.build("PyPlot")
using PyPlot

Pkg.add("ForwardDiff")
Pkg.add("HDF5")
Pkg.add("MAT")
Pkg.add("JLD")

using ForwardDiff
using MAT
using JLD

