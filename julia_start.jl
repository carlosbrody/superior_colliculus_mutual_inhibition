# Run this file in julia, and it will add and build all the necessary packages for cost function minimization

import Pkg

# === On Google VMs, we use the following Python install, but on Spock
# we're going to use Julia's PyCall's own Conda install of Python, so
# we comment out:
#
# ENV["PYTHON"]="/usr/bin/python"
#

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

using ForwardDiff
using MAT
using JLD
using MultivariateStats
using Clustering

# =====  the following for GPUs but commenting out for now
#
# if isfile("/usr/bin/lspci") && contains(readchomp(`lspci`), "NVIDIA")
#    Pkg.add("CUDAnative")
#    Pkg.add("CuArrays")
#    using CUDAnative
# end
