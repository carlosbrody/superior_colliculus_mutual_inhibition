using PyPlot
using ForwardDiff
using DiffResults
using MAT
using Random
using Optim
using JLD
using DelimitedFiles

JLD.translate("Core.Bool", "oldbool")

# G = load("farm_C30_Farms_C30_spock-brody01-03_0064.jld")

push!(LOAD_PATH, ".")

using Revise
using GradientUtils
using RateNetworks
using ProAnti
using TanhWalls

##

model_params[:opto_times] =["trial_start", "trial_start"]

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
