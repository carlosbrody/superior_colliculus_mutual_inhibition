

include("pro_anti_3node.jl")



mypars = Dict(
:init_add               =>          0,
:const_add              =>          0,
:noise                  =>          Any[],
:input                  =>          0,
:start_anti             =>          [-0.5, -0.5, -0.5],
:start_pro              =>          [-0.5, -0.5, -0.5],
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
:vW_2P                  =>          0,   # This and the following ones will be trained parameters
:vW_2A                  =>          0,   
:hW_P                   =>          0,
:sW_P                   =>          0,
:sW_A                   =>          0,
:sigma                  =>          0,
:constant_excitation    =>          0,
:target_period_excitation       =>  0,
:target_period_anti_excitation  =>  0,
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
                              "trial_start-0.1"     "trial_start-0.2"  
                              ],
:opto_targets     =>   [
                        0.9   0.7  
                        ],
)


search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW_2P  =>                   [mypars[:vW_2P],                    [-1,     1],  [-3,   3]], 
:vW_2A  =>                   [mypars[:vW_2A],                    [-1,     1],  [-3,   3]], 
:hW_P   =>                   [mypars[:hW_P],                     [-1,     1],  [-3,   3]],
:sW_P   =>                   [mypars[:sW_P],                     [0,      1],  [0,    3]],
:sW_A   =>                   [mypars[:sW_A],                     [0,      1],  [0,    3]],
:sigma  =>                   [mypars[:sigma],                    [0.05, 0.3],  [-2,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-2,     2],  [-30, 30]], 
:target_period_excitation => [mypars[:target_period_excitation], [-1,     1],  [-30  30]],
:target_period_anti_excitation   => [mypars[:target_period_anti_excitation],   [-1,     1],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-1,     1],  [-30, 30]],
#:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.5,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.5,],  [0,   30]],
)

#:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],

bbox = Dict()
for k in keys(search_conditions)
    bbox = merge(bbox, Dict(k=>search_conditions[k][3]))
end

