ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
using PyPlot

Pkg.add("ForwardDiff")
Pkg.add("HDF5")
Pkg.add("MAT")
