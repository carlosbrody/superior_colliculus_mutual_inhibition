## ====  testing OptimizationUtils

using Printf
using OptimizationUtils
using PyPlot


sr = Int64(round(time()*1000))
# sr = 1510162002784 # For these values of sr
# sr = 1509561656447 # when start_eta=1, the threshold quickly goes very positive and the minimization gets stuck
#
# sr = 1510164239381   # For this value, it gets stuck at a very small inverse slope param... don't know why
Random.seed!(sr)


npoints = 1000; # srand(400)
args = ["baseline", "amplitude", "threshold", "slope"]

# Generating values for our four params:
params = [1 5 0.5 0.8]

# Make some points and plot them
x = rand(npoints, 1)*6 .- 3
y = params[1] .+ params[2]*0.5*(tanh.((x.-params[3])./params[4]).+1) .+ randn(npoints,1)*2
figure(1); clf();
plot(x, y, ".")

# Starting values for the four params. Plot the corresponding curve they generate
seed = [8, 3.1, 0, 0.02]
xx = -3:0.01:3
plot(xx, seed[1] .+ seed[2]*0.5*(tanh.((xx.-seed[3])/seed[4]).+1), "g-")

# Cost function.  First output is the scalar
# that will be minimized, and we also returns a second output whose trajectory will be stashed
# by bbox in ftraj as a diagnostic during the minimization.
function myJJ(x, y; baseline=0, amplitude=1, threshold=0, slope=1, do_plot=false, fignum=1, clearfig=true)

    if do_plot
        figure(fignum);
        if clearfig; clf(); end;
        xx = -3:0.01:3; x2=zeros(get_eltype((baseline,amplitude,threshold,slope)), size(xx,1), size(xx,2))
        for i=1:length(xx); x2[i]=xx[i]; end; xx= x2

        plot(x, y, ".")
        plot(xx, baseline .+ amplitude*0.5*(tanh.((xx.-threshold)/slope).+1), "r-")
    end

    yhat =  baseline .+ amplitude*0.5*(tanh.((x.-threshold)./slope).+1)
    err = yhat .- y
    return sum(err.*err), get_value(err)    # Note first output, the scalar to be minimized,
    # may be ForwardDiff Duals during the minimization, which is fine, so it can be differentiated.
    # The second one we use get_value to turn into regular Float64 array so it comes out readable.
end



if ~isdir("Trash"); mkdir("Trash"); end;  # we're going to put the iteration-step by iteration-step report file there

bbox = Dict(:baseline=>[-2, 10], :slope=>[0.001 5])
func = (;pars...) -> myJJ(x, y; do_plot=false, pars...)

stopping_func = (;cost=0, func_out=[], pars...) -> return cost<1500;   # Make that a high number and it'll stop early

opars, traj, cost, cpm_traj, ftraj = bbox_Hessian_keyword_minimization(seed, args, bbox, func,
    frac_cost_threshold = 0.5, stopping_function = stopping_func,
    verbose=true, verbose_level=1, verbose_file="tester.txt",
    softbox=true, start_eta=0.1, report_file="Trash/example_report.jld")

# Note that the gradient at step i of the minimization will be available as ftraj[1,i], the hessian will be
# in ftraj[2,i], and the error vector, which is the first of the extra outputs of myJJ(), will be in ftraj[3,i][1].
# In our example myJJ() produced only one extra output; a second extra output would be in ftraj[3,i][2], and so on.

# Plot the resulting curve, and report both final and generating params
figure(1);
plot(xx, opars[1] .+ opars[2]*0.5*(tanh.((xx.-opars[3])/opars[4]).+1), "r-")
[opars' ; params]
xlabel("x"); ylabel("y"); title("green is sigmoid with starting params, red is end")


figure(2); clf();
ax1 = subplot(2,1,1)
plot(cpm_traj[4,:], ".-")
plot(traj[2,2:end] - traj[2,1:end-1], ".-")
grid("on")
gca().set_xticklabels("")
legend(["expected cost change", "actual cost change"])

subplot(2,1,2)
plot((traj[2,2:end] - traj[2,1:end-1])./cpm_traj[4,1:end-1], ".-")
plot(traj[1,2:end]./traj[1,1:end-1], ".-")
grid("on")

legend(["actual/expected cost change", "fractional change in eta"])

println("Final costs were: ", traj[2,end-3:end]);



##   ====  testing using callbacks with package Optim

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
