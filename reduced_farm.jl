# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


include("pro_anti.jl")


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.



mypars = Dict(
:init_add               =>          0,
:const_add              =>          0,
:noise                  =>          Any[],
:input                  =>          0,
:seedrand               =>          1510340076445,
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
:nPro => 100,
:nAnti => 100,
)



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.



extra_pars = Dict(
:plot_list        =>   [], 
:plot_conditions  =>   true,
:verbose          =>   true,
:opto_periods     =>   String["trial_start" "trial_start-0.1"],
:opto_targets     =>   [0.9,  0.7],
:opto_times       =>   ["trial_start", "trial_start-0.1"],       # This one is for run_ntrials
:seedrand         =>   Int64(round(time()*1000)),  # 1510426840370  Works wonders with search_range = 0.01
:cbeta            =>   0.001,
:search_range     =>   0.8,
)


search_conditions = Dict(   # :param    default_start   search_box  bound_box
:vW     =>                   [mypars[:vW],                       [-0.5, 0.5],  [-3,   3]], 
:hW     =>                   [mypars[:hW],                       [-0.5, 0.5],  [-3,   3]],
:dW     =>                   [mypars[:dW],                       [-0.5, 0.5],  [-3,   3]],
:sW     =>                   [mypars[:sW],                       [0,    0.5],  [0,    3]],
:sigma  =>                   [0.11,                              [0.1,  0.2],  [0.1, 0.4]],
:constant_excitation      => [mypars[:constant_excitation],      [-1,     1],  [-2,   2]], 
:target_period_excitation => [mypars[:target_period_excitation], [0.001,0.5],  [0     4]],
:right_light_excitation   => [mypars[:right_light_excitation],   [0.05, 0.5],  [0.05, 4]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-0.5, 0.5],  [-2,   2]],
# :opto_strength            => [mypars[:opto_strength],            [0.7, 0.99],  [0,    1]],
)


search_range = extra_pars[:search_range]; 

fbasename = "FarmFields/farm_C2_"

while true
    args = []; seed = []; bbox = Dict()
    for k in keys(search_conditions)
        search_box = search_conditions[k][2]
        args = [args; String(k)]
        myseed = search_conditions[k][1] + search_range*(rand()-0.5)*diff(search_box); myseed = myseed[1]
        if myseed > search_box[2]; myseed = search_box[2]; end
        if myseed < search_box[1]; myseed = search_box[1]; end
        seed = [seed ;  myseed]
        # seed = [seed ; search_conditions[k][1]]
        bbox = merge(bbox, Dict(k => Array{Float64}(search_conditions[k][3])))
    end
    args = Array{String}(args)
    seed = Array{Float64}(seed)



    func =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, 
    merge(merge(mypars, extra_pars), Dict(params))...)[1]

    try
        pars, traj, cost, cpm_traj, ftraj = bbox_Hessian_keyword_minimization(seed, args, bbox, func, 
            start_eta = 0.1, tol=1e-12, 
            verbose=true, verbose_every=1, maxiter=1000)


        cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(10000, 10000; verbose=false, 
        make_dict(args, pars, merge(merge(mypars, extra_pars)))...)

        myfilename = next_file(fbasename, 4)
        myfilename = myfilename*".jld"
        # write file
        save(myfilename, Dict("nPro"=>mypars[:nPro], "nAnti"=>mypars[:nAnti], 
        "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, "bbox"=>bbox, 
        "pars"=>pars, "traj"=>traj, "cost"=>cost, "cpm_traj"=>cpm_traj, "ftraj"=>ftraj,
        "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))
    catch
        @printf("\n\nWhoopsety, unkown error!   Trying new random seed.\n\n")
    end
end


