# Run this file in julia, and it will add and build all the necessary packages for cost function minimization

ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
using PyPlot

Pkg.add("ForwardDiff")
Pkg.add("HDF5")
Pkg.add("MAT")
