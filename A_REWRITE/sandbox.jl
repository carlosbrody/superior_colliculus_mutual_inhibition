##  === Pulling in Reports from NEW-style training and parsing them

using DelimitedFiles
include("sixNodeSetup_C32.jl")
args = argsSeedBounder()[1]   # get the string list of model params



pullReports=false; if pullReports   # if there are new runs on the VMs, get them
    remotes = "spock"

    if remotes == "googleVMs"
        machines = ["proanti006", "proanti007", "proanti008",
                    "proanti009", "proanti010", "proanti011"]

        for i in machines
            run(`./pull_reports.sh $i`)
        end
    elseif remotes == "spock"
        run(`./pull_reports.sh`)
    else
        error("Don't know remotes $remotes")
    end

    # At the of this block, report files will have been put in ../../Reports
end
##
global costs  = fill(0.0, 0)   # the cost for each solution
global fnames = fill("",  0)   # the filename for each solution
global pvals  = fill(0.0, 0, length(args))  # model params for each solution
global srands = fill(0, 0)     # random seeds for each solution

# Now we'll parse each Report file to find where it ended up at.
u = readdir("../../Reports")
for i in filter(x -> startswith(x, "r6_spockRun3"), u)
    # Read the file into string s:
    fp = open("../../Reports/$i", "r")
    s = read(fp, String)
    close(fp)
    println(i)
    # Find the last cost report:
    z = findlast("OVERALL", s)
    if z != nothing  && findlast("random seed", s)!=nothing
        z2 = findfirst("cost=", s[z[end]:end])
        z3 = findfirst(",", s[z[end]:end])
        mycost = parse(Float64, s[z[end]+z2[end]:z[end]+z3[1]-2])
        # println("$i: mycost=$mycost")
        global costs  = vcat(costs, mycost)
        global fnames = vcat(fnames, i)

        z = findlast("ps=[", s)[end].+1
        n = findfirst("]", s[z+1:end])[1] + z-1
        global pvals = vcat(pvals, readdlm(Vector{UInt8}(s[z:n]), ','))

        sz1 = findlast("random seed", s)
        sz2 = findfirst('\n', s[sz1[end]:end])
        global srands = vcat(srands, parse(Int64, s[sz1[end]+2:sz1[end]+sz2-2]))
    end
    # println("Checked $i")
end

println("$(length(costs)) networks went into further training")

nguys = 35   # number of solutions with the best costs that we'll look at
colids = [5, 7, 21, 10, 22, 4, 18, 3] # args[colids] will be the weights collected into "merges"

# we only want to see the best cost solutions
p = sortperm(costs, rev=true)
costs = costs[p]
fnames = fnames[p]
pvals = pvals[p,:]
srands = srands[p]
display([fnames costs])

# merges will hold the chosen weights
merges = [reshape(args[colids],1,length(colids)) ; pvals[end-nguys:end, colids]]
merges = hcat(["cost" ; costs[end-nguys:end]], merges)
display(merges)

cost_threshold = -0.0001
u = findall(costs .<= cost_threshold)
println("\nFound $(length(u)) networks below cost threshold of $cost_threshold")
nguys = length(u)-1

matwrite("Solutions/solutions6.mat", Dict(
    "costs"=>costs[end-nguys:end],
    "fnames"=>fnames[end-nguys:end],
    "argnames"=>args,
    "paramvals"=>pvals[end-nguys:end,:],
    "randseeds"=>srands[end-nguys:end],
    "weightMatrixMap"=>makeWeightMatrix(mypars; asString=true)
    ))

# # After running the cell above, run this one to produce parameter scatterplot

pairs = [
    "vW_P1A1" "vW_P1A2"
    "dW_P1A1" "dW_P1A2";
    ]

coords = fill(0.0, length(u), 2)

for i=1:size(pairs,1)
    u1 = findfirst(args.==pairs[i,1])
    u2 = findfirst(args.==pairs[i,2])

    coords[:,i] = 0.5*(pvals[u,u1] .+ pvals[u,u2])
end

figure(11); clf()
markersize=12;
fontname="Helvetica"
fontsize=14
axlim=2.5

subplot(1,2,1)
bins = -2.2:0.2:2.2
hist(pvals[u, findfirst(args.=="hW_P1P1")])
# axis("square");
xlim(bins[1]-0.1, bins[end]+0.1)
vlines([0], ylim()[1], ylim()[2], linestyle="--", linewidth=1.5, color="k")
gca().tick_params(labelsize=fontsize)
xlabel("hW_P1P1", fontname=fontname, fontsize=fontsize);
title("Pro<-->Pro weight,\n$(length(u)) networks with cost < -0.0001", fontname=fontname, fontsize=fontsize)
axisMove(-0.025, 0)

subplot(1,2,2)
plot(coords[:,1], coords[:,2], "k.", markersize=12);
xlabel("mean vW_PA", fontname=fontname, fontsize=fontsize);
ylabel("mean dW_PA", fontname=fontname, fontsize=fontsize)
axis("square")
xlim([-axlim,axlim]); ylim([-axlim,axlim]);
vlines([0], ylim()[1], ylim()[2], linestyle="--", linewidth=1.5, color="k")
hlines([0], xlim()[1], xlim()[2], linestyle="--", linewidth=1.5, color="k")
plot([-axlim,axlim], [-axlim,axlim], "r--")
gca().tick_params(labelsize=fontsize)
title("$(length(u)) networks with cost < -0.0001", fontname=fontname, fontsize=fontsize)
axisMove(0.025, 0)

savefig2jpg("Plots/parameterScatterplot")

## After running the cell above, run this one to produce parameter scatterplot

mycolids = vcat(colids, [10,22])
collist = unique(vcat(mycolids, setdiff(1:25, mycolids)))

figure(10); clf()
for i=1:25
    nrows = 5; ncols =5;
    subplot(nrows, ncols, i)
    hist(pvals[end-nguys:end, collist[i]])
    myrow = Int64(floor((i-1)/nrows))+1
    axisMove(0, 0.025*(nrows/2+0.5 - myrow))
    axisWidthChange(0.9, lock="c"); axisHeightChange(0.9, lock="c")
    title(args[collist[i]])
end
savefig2jpg("Plots/parameterHistograms")

figure(11); clf();
pairs = [
    "dW_P1A1" "dW_P1A2";
    "vW_P1A1" "vW_P1A2"
    ]

for i=1:size(pairs,1)
    u1 = findfirst(args.==pairs[i,1])
    u2 = findfirst(args.==pairs[i,2])
    subplot(1,2,i)
    hist(0.5*(pvals[end-nguys:end, u1] .+ pvals[end-nguys:end, u2]))
    title("mean of $(args[u1]) and $(args[u2])")
end
savefig2jpg("Plots/parameterMeanHistograms")




## Running an old solution.

# Solutions 19, 21, and 23 are interesting
id = size(pvals,1)-24  # We choose a number between 1 and size(pvals,1) to identify
# the solution in pvals that we'll run with
opto_condition = 1   # 1 is control, 2 is rule_and_delay, 3 is target_period

params = pvals[id,:]
include("sixNodeSetup_C32.jl")
extra_pars[:seedrand] = srands[id]
extra_pars[:few_trials]                = 50
extra_pars[:many_trials]               = 1600

ntrials = 40
extra_pars[:nPro] = ntrials
extra_pars[:nAnti] = ntrials
# mypars[:opto_units] = 1:6

# params[15] = 0.01

# This will now run the same way as when the cost during training was calculated,
# meaning that costs[id] and answers["cost"] will be the same:
answers = JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=true, asDict=true,
    model_details=true, merge(make_dict(args, params, mypars), extra_pars)...)
display(answers)

# Plot the activity of some nodes
mid = id - size(pvals,1)+size(merges,1)   # id in the mergers
display(merges[[1;mid],1:end])
t = (1:length(answers["antiValls"][1,1][1,:,1]))*mypars[:dt] .- mypars[:dt]


figure(1); clf()
nodes2plot = [3 5 ; 2 4 ; 1 6]
for row=1:size(nodes2plot,1)
    for col=1:size(nodes2plot,2)
        node = nodes2plot[row, col]
        rule_and_delay_period_id = 1
        target_period_id = 1
        nrdp = length(mypars[:rule_and_delay_periods])
        ntp  = length(mypars[:target_periods])
        pid   = (rule_and_delay_period_id-1)*ntp + target_period_id
        proValls  = answers["proValls" ][opto_condition,pid]
        antiValls = answers["antiValls"][opto_condition,pid]

        subplot(size(nodes2plot,1),size(nodes2plot,2), (row-1)*size(nodes2plot,2)+col); # axisMove(0, 0.025*(2-node))
        plot(t, proValls[node,:,1:ntrials],  color="b", linewidth=0.5)
        plot(t, antiValls[node,:,1:ntrials], color="r", linewidth=0.5);
        if row<size(nodes2plot,1); remove_xtick_labels(); end

        radp = mypars[:rule_and_delay_periods][rule_and_delay_period_id]
        tagp = mypars[:target_periods][target_period_id]
        vlines([radp, radp+tagp], ylim()[1], ylim()[2], linewidth=0.5)
        if in(node, mypars[:ProNodeID])
            flavor = "Pro"
            title("$flavor node (id=$node), blue is Pro trials, red is Anti trials")
        else
            flavor = "Anti"
            title("$flavor node $(node-1), blue is Pro trials, red is Anti trials")
        end
    end
end

subplot(size(nodes2plot,1),size(nodes2plot,2),1)
vW_P1A1 = findfirst(args.=="vW_P1A1")
vW_P1A2 = findfirst(args.=="vW_P1A2")
t = text(1.5, 1.3, "Solution #$(size(pvals,1)-id),\nA1->P1=$(params[vW_P1A1])"*
    ",  A2->P1=$(params[vW_P1A2])",
    horizontalalignment="center", fontSize=15)

savefig2jpg("Plots/solution$(size(pvals,1)-id)")

display([args pvals[id,:]])

## === Linearizing around specific solutions/conditions/timepoints


id = size(pvals,1)-9  # We choose a number between 1 and size(pvals,1) to identify
# the solution in pvals that we'll run with

params = pvals[id,:]
include("sixNodeSetup_C32.jl")
extra_pars[:seedrand] = srands[id]

nPro  = 10;
nAnti = 10;
mypars[:rule_and_delay_period] = mean(mypars[:rule_and_delay_periods])
mypars[:target_period]         = mean(mypars[:target_periods])

my_params = merge(mypars, make_dict(args, params))

proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input =
    run_ntrials(nPro, nAnti; my_params...)




#

# my_params = make_dict(["rule_and_delay_period","target_period","post_target_period"],   [i,j,k])
# my_params = make_dict(["opto_times"], [reshape(opto_periods[nopto,:], 1, 2)], my_params)
# my_params = merge(Dict(model_params), my_params)  # my_params takes precedence


##  === Pulling in Reports from old-style training and parsing them

pullReports=false; if pullReports
    machines = ["proanti002", "proanti004", "proanti005"]

    for i in machines
        run(`./pull_reports.sh $i`)
    end
end
##
global costs  = fill(0.0, 0)
global fnames = fill("",  0)
global pvals  = fill(0.0, 0, 16)

u = readdir("../../Reports")
for i in u
    fp = open("../../Reports/$i", "r")
    s = read(fp, String)
    close(fp)

    z = findlast("OVERALL", s)
    if z != nothing
        z2 = findfirst("cost=", s[z[end]:end])
        z3 = findfirst(",", s[z[end]:end])
        mycost = parse(Float64, s[z[end]+z2[end]:z[end]+z3[1]-2])
        # println("$i: mycost=$mycost")
        global costs  = vcat(costs, mycost)
        global fnames = vcat(fnames, i)

        z = findlast("ps=[", s)[end].+1
        n = findfirst("]", s[z+1:end])[1] + z-1
        global pvals = vcat(pvals, readdlm(Vector{UInt8}(s[z:n]), ','))
    end
    # println("Checked $i")
end


nguys = 15
colids = [11, 16, 10]

p = sortperm(costs, rev=true)
costs = costs[p]
fnames = fnames[p]
pvals = pvals[p,:]
display([fnames costs])

merges = [reshape(args[colids],1,length(colids)) ; pvals[end-nguys:end, colids]]
merges = hcat(["cost" ; costs[end-nguys:end]], merges)
display(merges)
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
