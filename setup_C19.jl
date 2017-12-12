

include("pro_anti.jl")



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
:vW                     =>          0,   # This and the following ones will be trained parameters
:hW                     =>          0,
:dW                     =>          0,
:sW                     =>          0,
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
                               "target_start+0.016"  "target_end" ; 
                              ],
:opto_targets     =>   [
                        0.9   0.7  ;
                        0.85  0.5  ;
                        0.9   0.7  ;
                        ],
)


search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW     =>                   [mypars[:vW],                       [-0.5, 0.5],  [-3,   3]], 
:hW     =>                   [mypars[:hW],                       [-0.5, 0.5],  [-3,   3]],
:dW     =>                   [mypars[:dW],                       [-0.5, 0.5],  [-3,   3]],
:sW     =>                   [mypars[:sW],                       [0,    0.5],  [0,    3]],
:sigma  =>                   [mypars[:sigma],                    [0.1,  0.2],  [-2,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-1,     1],  [-30, 30]], 
:target_period_excitation => [mypars[:target_period_excitation], [0.001,0.5],  [-30  30]],
:right_light_excitation   => [mypars[:right_light_excitation],   [0.05, 0.5],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-0.5, 0.5],  [-30, 30]],
:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],
:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.2,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.2,],  [0,   30]],
)



bbox = Dict()
for k in keys(search_conditions)
    bbox = merge(bbox, Dict(k=>search_conditions[k][3]))
end

