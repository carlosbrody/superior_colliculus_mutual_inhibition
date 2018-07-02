
try
    println("Trying to add MAT package")
    Pkg.add("MAT")
    using MAT
    println("Success! MAT package added")
catch x
    println("FAILURE!!! MAT package")
    println(x)
end

try
    println("Trying to add DiffBase package")
    Pkg.add("DiffBase")
    using DiffBase
    println("Success! DiffBase package added")
catch x
    println("FAILURE!!! DiffBase package")
    println(x)
end

try
    println("Trying to add ForwardDiff package")
    Pkg.add("ForwardDiff")
    using ForwardDiff
    println("Success! ForwardDiff package added")
catch x
    println("FAILURE!!! ForwardDiff package")
    println(x)
end

try
    println("Trying to add PyCall package")
    Pkg.add("PyCall")
#    ENV["PYTHON"]=""
#    Pkg.build("PyCall")
    using PyCall
#    println("Trying to add PyPlot package")   
#    Pkg.add("PyPlot")
#    using PyPlot
    println("Success! PyCall and PyPlot package added")
catch x
    println("FAILURE!!! PyCall and PyPlot packages")
    println(x)
end

try
    println("Trying to add HDF5 package")
    Pkg.add("HDF5")
    using HDF5
    println("Success! HDF5 package added")
catch x
    println("FAILURE!!! HDF5 package")
    println(x)
end

try
    println("Trying to add JLD package")
    Pkg.add("JLD")
    using JLD
    println("Success! JLD package added")
catch x
    println("FAILURE!!! JLD package")
    println(x)
end

try
    println("Trying to add MultivariateStats package")
    Pkg.add("MultivariateStats")
    using MultivariateStats
    println("Success! MultivariateStats package added")
catch x
    println("FAILURE!!! MultivariateStats package")
    println(x)
end

try
    println("Trying to add Clustering package")
    Pkg.add("Clustering")
    using Clustering
    println("Success! Clustering package added")
catch x
    println("FAILURE!!! Clustering package")
    println(x)
end









Pkg.update()
