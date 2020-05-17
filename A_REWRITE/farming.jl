if length(ARGS)>0  &&  tryparse(Int64, ARGS[1]) != nothing
   my_run_number = parse(Int64, ARGS[1]);
else
   my_run_number = 1; # I am process my_run_number
end
if length(ARGS)>0  &&  tryparse(Int64, ARGS[2]) != nothing
   tot_n_runs    = parse(Int64, ARGS[2]);
else
   tot_n_runs = 1;   # I'm being run as part of tot_n_run processes
end

include("startup.jl")

mypars = Dict(
:init_add               =>          0,
:const_add              =>          0,
:noise                  =>          Any[],
:input                  =>          0,
:start_anti             =>          [-0.5, -0.5, -0.5, -0.5],
:start_pro              =>          [-0.5, -0.5, -0.5, -0.5],
:rule_and_delay_period  =>          1.2,
:rule_and_delay_periods =>          [1.2],
:target_period          =>          0.3,
:target_periods         =>          [0.3],
:post_target_period     =>          0.3,
:post_target_periods    =>          [0.3],
:right_light_pro_extra  =>          0,
:U_rest                 =>          0,
:theta                  =>          0.05,
:beta                   =>          0.5,
:g_leak                 =>          1,
:nsteps                 =>          301,
:dt                     =>          0.024,
:tau                    =>          0.09,
:symmetrized_W          =>          false,
:vW_PA                  =>          0,   # This and the following ones will be trained parameters
:vW_AP                  =>          0,
:hW_P                   =>          0,
:hW_A                   =>          0,
:dW_PA                  =>          0,
:dW_AP                  =>          0,
:sW_P                   =>          0,
:sW_A                   =>          0,
:sigma                  =>          0,
:constant_excitation    =>          0,
:target_period_excitation =>        0,
:right_light_excitation =>          0,
:const_pro_bias         =>          0,
:opto_strength          =>          0,
:pro_rule_strength      =>          0,
:anti_rule_strength     =>          0,
)


vars = [
    "vW_PA"                         0.1;
    "vW_AP"                         0.1;
    "hW_P"                          0.1;
    "hW_A"                          0.1;
    "dW_PA"                         0.1;
    "dW_AP"                         0.1;
    "sW_P"                          0.1;
    "sW_A"                          0.1;
    "sigma"                         0.05;
    "constant_excitation"           0.1;
    "target_period_excitation"      0.1;
    "right_light_excitation"        0.1;
    "const_pro_bias"                0.1;
    "opto_strength"                 0.1;
    "pro_rule_strength"             0.1;
    "anti_rule_strength"            0.1
    ]

args  = Array{String}(vars[:,1])
seed  = Array{Float64}(vars[:,2])

extra_pars = Dict(
:theta2           =>   0.15,
:theta1           =>   0.05,
:cbeta            =>   0.001,
:plot_list        =>   [],
:plot_conditions  =>   true,
:opto_times       =>   ["trial_start", "trial_start-0.1"],       # This one is a default for run_ntrials
:opto_periods     =>   String[
                              "trial_start-0.1"     "trial_start-0.2" ;
                              "target_start-0.4"    "target_start" ;
                               "target_start+0.016"  "target_end" ;
                              ],
:opto_targets     =>   [
                        0.9   0.7  ;
                        0.85  0.5  ;
                        0.9   0.7  ;
                        ],
)

extra_pars[:seedrand]                  = Int64(my_run_number*round(time()*1000))
extra_pars[:maxiter]                   = 1000
extra_pars[:testruns]                  = 10000
extra_pars[:few_trials]                = 50
extra_pars[:many_trials]               = 1600
extra_pars[:binarized_delta_threshold] = 0.1    # average frac correct must be within this of target
extra_pars[:anti_perf_delta]           = 0.05   # delay anti must be at least this worse off than control or choice anti
extra_pars[:pro_better_than_anti]      = true   # if true, in each condition pro hits must be > anti hits

extra_pars[:nPro]  = 50
extra_pars[:nAnti] = 50

using Random
Random.seed!(extra_pars[:seedrand])

##

G = load("farm_C30_Farms_C30_spock-brody01-03_0064.jld")

func =  x -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
    make_dict(G["args"], x, merge(mypars, extra_pars))...)[1]

g = x->ForwardDiff.gradient(func, x)
h = x->ForwardDiff.hessian(func, x)


## ======================================================
#
#   BOUNDING BOX
#
#  ======================================================

search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW_PA  =>                   [mypars[:vW_PA],                    [-1,     1],  [-3,   3]],
:vW_AP  =>                   [mypars[:vW_AP],                    [-1,     1],  [-3,   3]],
:hW_P   =>                   [mypars[:hW_P],                     [-1,     1],  [-3,   3]],
:hW_A   =>                   [mypars[:hW_A],                     [-1,     1],  [-3,   3]],
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

args   = Array{String}(undef, 0)
seed   = Array{Float64}(undef, 0)
lower  = Array{Float64}(undef, 0)
upper  = Array{Float64}(undef, 0)

for k in keys(search_conditions)
    sbox = search_conditions[k][2]
    bbox = search_conditions[k][3]

    global args  = vcat(args, string(k))
    global seed  = vcat(seed, sbox[1] + rand()*(diff(sbox)[1]))
    global lower = vcat(lower, bbox[1])
    global upper = vcat(upper, bbox[2])
end

ResultsStash = ["seedrand" "cost" "minimizer" ;
                 0        func(seed)  [seed]]



result = optimize(
    func, g, # h, seems overall faster without?
    lower, upper,
    seed, # NewtonTrustRegion();  seems overall faster without?
    Fminbox(),
    Optim.Options(store_trace=true, show_trace=true, time_limit=30000);
    inplace=false);

println("With seedrand ", extra_pars[:seedrand], " true cost was ",
    func(Optim.minimizer(result)))

ResultsStash = vcat(ResultsStash,
    [extra_pars[:seedrand] Optim.minimum(result) [Optim.minimizer(result)]])

println("\nResultsStash:")
display(ResultsStash)
println()

println("\nJJ:")
display(JJ(extra_pars[:nPro], extra_pars[:nAnti]; asDict=true, verbose=false,
    make_dict(G["args"], Optim.minimizer(result), merge(mypars, extra_pars))...))

println()
