

using ProAnti

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


mypars = Dict(
:init_add               =>          0,
:const_add              =>          0,
:noise                  =>          Any[],
:input                  =>          0,
:start_anti             =>          [-0.5, -0.5, -0.5, -0.5],
:start_pro              =>          [-0.5, -0.5, -0.5, -0.5],
:rule_and_delay_period  =>          1.2,
:rule_and_delay_periods =>          [1.0 1.2],
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
:sigma                  =>          0,
:constant_excitation    =>          0,
:target_period_excitation =>        0,
:right_light_excitation =>          0,
:const_pro_bias         =>          0,
:opto_strength          =>          0,
:pro_rule_strength      =>          0,
:anti_rule_strength     =>          0,
:AntiNodeID             =>          [2,3,4,5],
:ProNodeID              =>          [1,6],
:RightNodeID            =>          [1,2,3],
:LeftNodeID             =>          [4,5,6],
:start_anti             =>          -0.5*ones(6),
:start_pro              =>          -0.5*ones(6),
:opto_units             =>          1:6
)

for i in W6
    mypars[Symbol(i)] = 0
end


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


search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW_P1A1  =>                 [mypars[:vW_P1A1],                  [-1,     1],  [-3,   3]],
:vW_P1A2  =>                 [mypars[:vW_P1A2],                  [-1,     1],  [-3,   3]],
:vW_A1P1  =>                 [mypars[:vW_A1P1],                  [-1,     1],  [-3,   3]],
:vW_A2P1  =>                 [mypars[:vW_A2P1],                  [-1,     1],  [-3,   3]],
:hW_P1P1  =>                 [mypars[:hW_P1P1],                  [-1,     1],  [-3,   3]],
:cW_A2A1  =>                 [mypars[:cW_A2A1],                  [-1,     1],  [-3,   3]],
:cW_A1A2  =>                 [mypars[:cW_A1A2],                  [-1,     1],  [-3,   3]],
:hW_A2A1  =>                 [mypars[:hW_A2A1],                  [-1,     1],  [-3,   3]],
:hW_A1A2  =>                 [mypars[:hW_A1A2],                  [-1,     1],  [-3,   3]],
:hW_A2A2  =>                 [mypars[:hW_A2A2],                  [-1,     1],  [-3,   3]],
:hW_A1A1  =>                 [mypars[:hW_A1A1],                  [-1,     1],  [-3,   3]],
:dW_P1A1  =>                 [mypars[:dW_P1A1],                  [-1,     1],  [-3,   3]],
:dW_P1A2  =>                 [mypars[:dW_P1A2],                  [-1,     1],  [-3,   3]],
:dW_A1P1  =>                 [mypars[:dW_A1P1],                  [-1,     1],  [-3,   3]],
:dW_A2P1  =>                 [mypars[:dW_A2P1],                  [-1,     1],  [-3,   3]],
:sW_P1    =>                 [mypars[:sW_P1],                    [0,      1],  [0,    3]],
:sW_A1    =>                 [mypars[:sW_A1],                    [0,      1],  [0,    3]],
:sW_A2    =>                 [mypars[:sW_A2],                    [0,      1],  [0,    3]],
:sigma  =>                   [mypars[:sigma],                    [0.05, 0.3],  [-2,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-2,     2],  [-30, 30]],
:target_period_excitation => [mypars[:target_period_excitation], [-1,     1],  [-30  30]],
:right_light_excitation   => [mypars[:right_light_excitation],   [-1,     1],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-1,     1],  [-30, 30]],
:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],
:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.5,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.5,],  [0,   30]],
)



global bbox = Dict()
for k in keys(search_conditions)
    global bbox = merge(bbox, Dict(k=>search_conditions[k][3]))
end

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
