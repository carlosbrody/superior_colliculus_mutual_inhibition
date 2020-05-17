# Run this file in julia, and it will add and build all the necessary packages for cost function minimization

ENV["PYTHON"]="/usr/bin/python"

import Pkg

Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
Pkg.build("PyPlot")
using PyPlot

Pkg.add("ForwardDiff")
Pkg.add("DiffResults")

Pkg.add("HDF5")
Pkg.add("MAT")
Pkg.add("JLD")
Pkg.add("MultivariateStats")
Pkg.add("Clustering")
Pkg.add("Optim")
Pkg.add("Revise")

using ForwardDiff
using MAT
using JLD
using MultivariateStats
using Clustering
using Optim
using Revise

if isfile("/usr/bin/lspci") && contains(readchomp(`lspci`), "NVIDIA")
    Pkg.add("CUDAnative")
    Pkg.add("CuArrays")
    using CUDAnative
end
