# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.



include("pro_anti.jl")


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


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
:anti_rule_strength     =>          0.054,
:U_rest                 =>          0,
:theta                  =>          0.05,
:beta                   =>          0.5,
:g_leak                 =>          1,
:nsteps                 =>          301,
:dt                     =>          0.024,
:tau                    =>          0.09,
:right_light_excitation =>          0.49924152955481954,
:opto_strength          =>          0.85,
:opto_periods           =>          String[
                                    "trial_start" "trial_start"; 
                                    "trial_start" "trial_end"; 
                                    "trial_start" "target_start/2"; 
                                    "target_start/2" "target_start"; 
                                    "target_start" "target_end"],
:opto_targets          =>           [
                                    0.9 0.7; 
                                    0.9 0.5; 
                                    0.9 0.7; 
                                    0.9 0.5; 
                                    0.9 0.7],
:theta2 => 0.15,
:theta1 => 0.05,
:sigma => 0.01,
:cbeta => 0.04,
:sW => 0.6416875048295452,
:hW => 0.054701836208134846,
:dW => 0.1267124266934907,
:vW => -1.588850577499782,
:constant_excitation => -0.37242520737694207,
:const_pro_bias => 0.04366897857834884,
:target_period_excitation => 0.15315254453690974,
:right_light_pro_extra => 0,
:pro_rule_strength => 0.05,
:nPro => 200,
:nAnti => 200,
)


extra_pars = Dict(
:plot_list        =>   [], 
:plot_conditions  =>   true,
:verbose          =>   true,
:opto_periods     =>   String["trial_start" "trial_start-0.1"],
:opto_targets     =>   [0.9,  0.7],
:opto_times       =>   ["trial_start", "trial_start-0.1"],       # This one is for run_ntrials
:cbeta            =>   0.001,
:search_range     =>   1.2,
)

search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW     =>                   [mypars[:vW],                       [-0.5, 0.5],  [-5,   5]], 
:hW     =>                   [mypars[:hW],                       [-0.5, 0.5],  [-5,   5]],
:dW     =>                   [mypars[:dW],                       [-0.5, 0.5],  [-5,   5]],
:sW     =>                   [mypars[:sW],                       [0,    0.5],  [-5,   5]],
:sigma  =>                   [0.11,                              [0.1,  0.2],  [-2,   2]],
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
bbox



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


README = """
Attempt to optimize the C17 farms, starting from all ending points that had
test costs <= 0.
"""

source_dir = "Available_C17_Farms"
optim_dir  = "Optimized_C17_Farms"
cost_threshold = 0   # only work with test costs less than this

if ~isdir(optim_dir); mkdir(optim_dir); end


f = readdir(source_dir)

while length(f) > 0
    mypars, extra_pars, args, seed, test_cost = load(source_dir * "/" * f[1], 
        "mypars", "extra_pars", "args", "pars3", "cost")
    rm(source_dir * "/" * f[1])

    if test_cost <= 0
        mypars[:nPro]  = 1600
        mypars[:nAnti] = 1600

        maxiter1 = 2 # 1000;   # for func1, the regular search
        testruns = 10000;  # Number of trials for evaluating the results of the model. 10000 is a good number 

        func1 =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, 
            merge(merge(mypars, extra_pars), Dict(params))...)[1]

        try
            pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(seed, 
                args, bbox, func1, 
                start_eta = 0.01, tol=1e-12, 
                verbose=true, verbose_every=1, maxiter=maxiter1)

            # evaluate the result with many trials, for accuracy
            cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(testruns, testruns; verbose=false, 
                make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)

            # Write out the results
            myfilename = optim_dir * "/" * f[1]

            @printf("\n\n ****** writing to file %s *******\n\n", myfilename)

            # write file
            save(myfilename, Dict("README"=>README, "nPro"=>mypars[:nPro], "nAnti"=>mypars[:nAnti], 
                "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, 
                "pars3"=>pars3, "traj3"=>traj3, "cost3"=>cost3, "cpm_traj3"=>cpm_traj3, "ftraj3"=>ftraj3,
                "cost"=>cost, "cost1s"=>cost1s, "cost2s"=>cost2s,
                "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))

        catch y
            # Interrupts should not get caught:
            if isa(y, InterruptException); throw(InterruptException()); end

            # Other errors get caught and a warning is issued but then we run again
            @printf("\n\nWhoopsety, unkown error!\n\n");
            @printf("Error was :\n"); print(y); @printf("\n\nTrying new random seed.\n\n")
        end
    end
    
    f = readdir(source_dir)
end


