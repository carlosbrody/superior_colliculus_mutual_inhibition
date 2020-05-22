##

using Optim

global latest_minimizer
global latest_minimum = Inf

println("\n\n")

function fu(x)
    answer = sum((x.-3).^2) + sum(x.^3)
    if answer < latest_minimum
        global latest_minimum = answer
        global latest_minimizer = x
    end
    return answer
end

g = x->ForwardDiff.gradient(fu, x)

function mycallback(x)
    if typeof(x) <: Array
        x = x[end]
    end
    println("$(x.iteration): pars=", get_value(latest_minimizer))
    # println("Cost now is ", x.value)
    # println("fu(latest_minimizer) = ", fu(latest_minimizer))
    return false
end

println("---")
result = optimize(
    fu, g, # seems overall faster without?
    randn(4), # NewtonTrustRegion(), # seems overall faster without?
    Optim.Options(show_trace=true, callback=mycallback, iterations=10);
    inplace=false)
##
#callback = mycallback,
#    store_trace=true, show_trace=true,
#    iterations=extra_pars[:secondPassNIter

include("commonSetup.jl")
args = argsSeedBounder()[1]

seed = Array{Float64}(readCostsFile(maxCost=-0.0002)[2,3:end])

fullpars = make_dict(args, seed, merge(model_params, extra_pars))

W6 = [
    "sW_P1"
    "sW_A1"
    "sW_A2"
    "vW_A1P1"
    "vW_A2P1"
    "dW_A1P1"
    "dW_A2P1"
    "hW_P1P1"
    "vW_P1A1"
    "cW_A2A1"
    "hW_A1A1"
    "hW_A2A1"
    "dW_P1A1"
    "vW_P1A2"
    "cW_A1A2"
    "hW_A1A2"
    "hW_A2A2"
    "dW_P1A2"
    ]


fullpars = merge(Dict(:AntiNodeID=>[2,3,4,5], :ProNodeID=>[1,6],
        :RightNodeID=>[1,2,3], :LeftNodeID=>[4,5,6]), fullpars)
fullpars = make_dict(W6, [rand(3); randn(15)], fullpars)
fullpars[:start_anti] = -0.5*ones(6)
fullpars[:start_pro]  = -0.5*ones(6)

makeWeightMatrix(fullpars; asString=false)

# run_ntrials(5, 5; fullpars...)


# ======================================================
#
#   BOUNDING BOX
#
#  ======================================================

search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW_P1A1  =>                 [mypars[:vW_P1A1],                  [-1,     1],  [-3,   3]],
:vW_P1A2  =>                 [mypars[:vW_P1A2],                  [-1,     1],  [-3,   3]],
:vW_A1P1  =>                 [mypars[:vW_A1P1],                  [-1,     1],  [-3,   3]],
:vW_A2P1  =>                 [mypars[:vW_A2P1],                  [-1,     1],  [-3,   3]],
:hW_P1    =>                 [mypars[:hW_P1],                    [-1,     1],  [-3,   3]],
:hW_A2A1  =>                 [mypars[:hW_A2A1],                  [-1,     1],  [-3,   3]],
:hW_A1A2  =>                 [mypars[:hW_A1A2],                  [-1,     1],  [-3,   3]],
:hW_A2   =>                  [mypars[:hW_A2],                    [-1,     1],  [-3,   3]],
:dW_PA  =>                   [mypars[:dW_PA],                    [-1,     1],  [-3,   3]],
:dW_AP  =>                   [mypars[:dW_AP],                    [-1,     1],  [-3,   3]],
:sW_P   =>                   [mypars[:sW_P],                     [0,      1],  [0,    3]],
:sW_A   =>                   [mypars[:sW_A],                     [0,      1],  [0,    3]],
:sigma  =>                   [mypars[:sigma],                    [0.05, 0.3],  [0,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-2,     2],  [-30, 30]],
:target_period_excitation => [mypars[:target_period_excitation], [-1,     1],  [-30  30]],
:right_light_excitation   => [mypars[:right_light_excitation],   [-1,     1],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-1,     1],  [-30, 30]],
:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],
:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.5,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.5,],  [0,   30]],
)

##


using Optim
using TanhWalls

function f(x)
   return sum(x.^2)
end

bfunc, old2new, new2old = boundedVersion(f, [1 3; -100 100])

g = x -> ForwardDiff.gradient(bfunc, x)

result = optimize(
    bfunc, g, # h, seems overall faster without?
    old2new([2.0,2]), # NewtonTrustRegion();  seems overall faster without?
    Optim.Options(store_trace=true, show_trace=true, time_limit=500);
    inplace=false);
