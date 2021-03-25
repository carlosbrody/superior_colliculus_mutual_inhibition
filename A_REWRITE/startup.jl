#
# Within the package manager, if using a new Julia:
#   add PyCall PyPlot ForwardDiff DiffResults MAT Random Optim JLD DelimitedFiles Dates Statistics ArgParse Debugger Revise
#   update Revise  # don't forget to update Revise! otherwise Julia terminates when using Revise

if !occursin("spock", String(chomp(read(`hostname`, String))))
    using PyPlot
    PyPlot.rc("font", family="Helvetica Neue")
end
using ForwardDiff
using DiffResults
using MAT
using Random
using Optim
using JLD
using DelimitedFiles
using Dates
using Statistics
using ArgParse

JLD.translate("Core.Bool", "oldbool")

println("\n     finished loading system modules in startup.jl ---")

# G = load("farm_C30_Farms_C30_spock-brody01-03_0064.jld")
##
if !in(".", LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using Revise
using GeneralUtils

## a

using GradientUtils
using RateNetworks
using ProAnti
using TanhWalls
using ResultsAnalysis
using ConstrainedParabolicMinimization
using OptimizationUtils

## a

model_params[:opto_times] =["trial_start", "trial_start"]

println("     finished loading user modules in startup.jl ---")
println("\n--- finished loading startup.jl ---\n")

# fun = x -> begin;
#    Random.seed!(20);
#    run_ntrials(10, 0; make_dict(["sW"], [x], model_params)...)[1][1];
# end
#
# ForwardDiff.gradient(fun, [0.4])
#
# ##
#
# fun = x -> begin;
#    Random.seed!(20);
#    JJ(50, 50; make_dict(["sW"], [x], model_params)...)[1];
# end
#
# ##
# ForwardDiff.gradient(fun, [0.4])
