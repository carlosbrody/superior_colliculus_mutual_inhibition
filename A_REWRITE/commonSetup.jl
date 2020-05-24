"""
This script sets my_run_number, tot_n_runs, and defines

    mypars
    extrapars
    search_conditions

    extra_pars[:seedrand] is not set in this script
    func() is not set in this script
"""


if length(ARGS)>0  &&  tryparse(Int64, ARGS[1]) != nothing
   my_run_number = parse(Int64, ARGS[1]);
else
   my_run_number = 1; # I am process my_run_number
end
if length(ARGS)>1  &&  tryparse(Int64, ARGS[2]) != nothing
   tot_n_runs    = parse(Int64, ARGS[2]);
else
   tot_n_runs = 1;   # I'm being run as part of tot_n_run processes
end


mypars = Dict(
:init_add               =>          0,
:const_add              =>          0,
:noise                  =>          Any[],
:input                  =>          0,
:start_anti             =>          [-0.5, -0.5, -0.5, -0.5],
:start_pro              =>          [-0.5, -0.5, -0.5, -0.5],
:rule_and_delay_period  =>          1.2,
:rule_and_delay_periods =>          [1.15 1.2],
:target_period          =>          0.6,
:target_periods         =>          [0.45 0.6],
:post_target_period     =>          0,
:post_target_periods    =>          [0],
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
                               "target_start+0.016"  "trial_end" ;
                              ],
:opto_targets     =>   [
                        0.9   0.7  ;
                        0.85  0.5  ;
                        0.9   0.7  ;
                        ],
)

##

# Enough trials, iters for a real run:
extra_pars[:nPasses]                   =  8       # of pass blocks below

extra_pars[:pass1NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass1NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass1CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass1RnD]                  = [1.2]    # rule and delay period range

extra_pars[:pass2NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass2NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass2CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass2RnD]                  = [1.15  1.2]    # rule and delay period range

extra_pars[:pass3NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass3NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass3CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass3RnD]                  = [1.1   1.2]    # rule and delay period range

extra_pars[:pass4NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass4NIter]                = 200      # maximum iterations in first pass
extra_pars[:pass4CostThreshold]        = 0        # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass4RnD]                  = [1.05   1.2]    # rule and delay period range

extra_pars[:pass5NTrials]              = 25       # number of trials to use in first pass
extra_pars[:pass5NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass5CostThreshold]        = -0.0001  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass5RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass6NTrials]              = 100      # number of trials to use in first pass
extra_pars[:pass6NIter]                = 400      # maximum iterations in first pass
extra_pars[:pass6CostThreshold]        = -0.0001  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass6RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass7NTrials]              = 400       # number of trials to use in first pass
extra_pars[:pass7NIter]                = 400       # maximum iterations in first pass
extra_pars[:pass7CostThreshold]        = -0.00015  # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass7RnD]                  = [1.0   1.2]    # rule and delay period range

extra_pars[:pass8NTrials]              = 1600       # number of trials to use in first pass
extra_pars[:pass8NIter]                = 400        # maximum iterations in first pass
extra_pars[:pass8CostThreshold]        = -0.00028   # maximum cost threshold for a first pass run to seed a second pass run
extra_pars[:pass8RnD]                  = [1.0   1.2]    # rule and delay period range

# Few trials, iters for testing:
# extra_pars[:nPasses]                   =  3       # of pass blocks below
#
# extra_pars[:pass1NTrials]              = 5       # number of trials to use in first pass
# extra_pars[:pass1NIter]                = 3      # maximum iterations in first pass
# extra_pars[:pass1CostThreshold]        = 2        # maximum cost threshold for a first pass run to seed a second pass run
# extra_pars[:pass1RnD]                  = [1.075 1.125]    # rule and delay period range
#
# extra_pars[:pass2NTrials]              = 5       # number of trials to use in first pass
# extra_pars[:pass2NIter]                = 3      # maximum iterations in first pass
# extra_pars[:pass2CostThreshold]        = 2        # maximum cost threshold for a first pass run to seed a second pass run
# extra_pars[:pass2RnD]                  = [1.05  1.1]    # rule and delay period range
#
# extra_pars[:pass3NTrials]              = 5       # number of trials to use in first pass
# extra_pars[:pass3NIter]                = 3      # maximum iterations in first pass
# extra_pars[:pass3CostThreshold]        = 2        # maximum cost threshold for a first pass run to seed a second pass run
# extra_pars[:pass3RnD]                  = [1.0   1.2]    # rule and delay period range




extra_pars[:binarized_delta_threshold] = 0.1    # average frac correct must be within this of target
extra_pars[:anti_perf_delta]           = 0.05   # delay anti must be at least this worse off than control or choice anti
extra_pars[:pro_better_than_anti]      = true   # if true, in each condition pro hits must be > anti hits
extra_pars[:maxiter]                   = 1000
extra_pars[:testruns]                  = 10000

extra_pars[:nPro]  = 50
extra_pars[:nAnti] = 50





# ======================================================
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
:sigma  =>                   [mypars[:sigma],                    [0.05, 0.3],  [-2,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-2,     2],  [-30, 30]],
:target_period_excitation => [mypars[:target_period_excitation], [-1,     1],  [-30  30]],
:right_light_excitation   => [mypars[:right_light_excitation],   [-1,     1],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-1,     1],  [-30, 30]],
:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],
:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.5,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.5,],  [0,   30]],
)

"""
   argsSeedBounder()

   This function uses the global search_conditions to define its return
   values.

   = RETURNS:

   - args     A vector of N strings defining the argumenrs

   - seed     A vector of N Float64s, each correspinding to its arg

   - bounder  An N-by-2 matrix, first column lower bounds, second column upper bounds.

"""
function argsSeedBounder()
   args    = Array{String}(undef, 0)
   seed    = Array{Float64}(undef, 0)
   bounder = Array{Float64}(undef, 0, 2)

   for k in keys(search_conditions)
      sbox = search_conditions[k][2]
      bbox = search_conditions[k][3]

      args    = vcat(args, string(k))
      seed    = vcat(seed, sbox[1] + rand()*(diff(sbox)[1]))
      bounder = vcat(bounder, bbox[:]')
   end

   return args, seed, bounder
end
