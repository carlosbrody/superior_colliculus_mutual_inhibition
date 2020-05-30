
println("Hey, im running an example script")
using MAT
using ForwardDiff
using DiffResults
using PyCall
using JLD
using MultivariateStats
using Clustering


myfilename = "test.mat"
test = 1

println("I loaded five packages into the environment, and made some local variables")

matwrite(myfilename, Dict("test"=>test))

println("I saved a file, trying reading with x=matread(\"test\")")


