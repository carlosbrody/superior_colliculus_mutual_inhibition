# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


README_TOP = """

Code for doing a minimization (no opto) in which we first do a quick search
with few trials, not on %correct but instead on specific differences between 
V[1] and V[4], for example for V[1]-V[4]=0.1 on Pro trials and =-0.1 on Anti 
trials.

We stop that quick pre-search when the binarized hits for Pro and for Anti are
both >= 0.7; and it is at those params that we then run the usual minimization.

This doesn't make 100% of minimization runs successful, but it does seem to
increase the hit rate. 

"""

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
:nPro => 100,
:nAnti => 100,
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
:vW     =>                   [mypars[:vW],                       [-0.5, 0.5],  [-3,   3]], 
:hW     =>                   [mypars[:hW],                       [-0.5, 0.5],  [-3,   3]],
:dW     =>                   [mypars[:dW],                       [-0.5, 0.5],  [-3,   3]],
:sW     =>                   [mypars[:sW],                       [0,    0.5],  [0,    3]],
:sigma  =>                   [0.11,                              [0.1,  0.2],  [-2,   2]],
:constant_excitation      => [mypars[:constant_excitation],      [-1,     1],  [-30, 30]], 
:target_period_excitation => [mypars[:target_period_excitation], [0.001,0.5],  [-30  30]],
:right_light_excitation   => [mypars[:right_light_excitation],   [0.05, 0.5],  [-30, 30]],
:const_pro_bias           => [mypars[:const_pro_bias],           [-0.5, 0.5],  [-30, 30]],
:opto_strength            => [mypars[:opto_strength],            [0,      1],  [0,    1]],
:pro_rule_strength        => [mypars[:pro_rule_strength],        [0,   0.2,],  [0,   30]],
:anti_rule_strength       => [mypars[:anti_rule_strength],       [0,   0.2,],  [0,   30]],
)




# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


extra_pars[:seedrand] = Int64(round(1000*time()))   
srand(extra_pars[:seedrand])
search_range = extra_pars[:search_range]; 

README = """

Farm C17: Like C16, with opto_strength and pro_rule_stregth and anti_rule_strength as 
trainable parameters, but with reasonable bounds on the appropriate parameters,
not no bounds at all.

"""

extra_pars[:opto_periods]  = String[
    "trial_start-0.1"     "trial_start-0.2" ; 
    "target_start-0.4"    "target_start" ; 
    "target_start+0.016"  "target_end" ; 
]

extra_pars[:opto_targets] = [
    0.9   0.7  ;
    0.85  0.5  ;
    0.9   0.7  ;
]

ftraj2 = []; cost2 = [];

farmdir = "../Farms" * chomp(readstring(`hostname`))[end-2:end]
if !isdir(farmdir); mkdir(farmdir); end
fbasename = farmdir * "/" * "farm_C17_"

reports_dir = "../Reports" * chomp(readstring(`hostname`))[end-2:end]
if !isdir(reports_dir); mkdir(reports_dir); end

my_run_number = tot_n_runs = 1
try
    if ~isnull(tryparse(Int64, ARGS[1])); my_run_number = parse(Int64, ARGS[1]); 
    else                                  my_run_number = 1; 
    end
    if ~isnull(tryparse(Int64, ARGS[2])); tot_n_runs    = parse(Int64, ARGS[2]); 
    else                                  tot_n_runs = 1; 
    end
catch
end

report_file = reports_dir * "/" * @sprintf("report_out_%d", my_run_number)

append_to_file(report_file, @sprintf("\n\n\nStarting with random seed %d\n\n\n", extra_pars[:seedrand]))

while true
    
    append_to_file(report_file, @sprintf("\n\n--- new run -- %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))
    args = []; seed = []; bbox = Dict()
    for k in keys(search_conditions)
        search_box = search_conditions[k][2]
        args = [args; String(k)]
        # --- search within the indicated search range:
        myseed = search_conditions[k][1] + search_range*(rand()-0.5)*diff(search_box); myseed = myseed[1]
        if myseed > search_box[2]; myseed = search_box[2]; end
        if myseed < search_box[1]; myseed = search_box[1]; end
        seed = [seed ;  myseed]
        # --- No search, just start at the indicated position:
        # seed = [seed ; search_conditions[k][1]]
        # --- search within the full indicated search box
        # seed = [seed ; rand()*diff(search_box) + search_box[1]]
        bbox = merge(bbox, Dict(k => Array{Float64}(search_conditions[k][3])))
    end
    args = Array{String}(args)
    seed = Array{Float64}(seed)


    maxiter1 = 2;   # for func1, the regular search
    testruns = 10;  # Number of trials for evaluating the results of the model. 10000 is a good number 


    # Make sure to keep the noise frozen over the search, meaning JJ() needs the seedrand parameter
    func1 =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, verbose_file=report_file, 
        merge(merge(mypars, extra_pars), Dict(params))...)[1]
    
    try
        pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(seed, 
            args, bbox, func1, 
            start_eta = 0.01, tol=1e-12, verbose_file=report_file,
            verbose=true, verbose_every=1, maxiter=maxiter1)
            
        # evaluate the result with many trials, for accuracy
        cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(testruns, testruns; verbose=false, 
            make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)
        
        # Write out the results
        myfilename = next_file(fbasename, 4)
        myfilename = myfilename*".jld"

        append_to_file(report_file, @sprintf("\n\n ****** writing to file %s *******\n\n", myfilename))
        
        # write file
        save(myfilename, Dict("README"=>README, "nPro"=>mypars[:nPro], "nAnti"=>mypars[:nAnti], 
            "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, "bbox"=>bbox, 
            "search_conditions"=>search_conditions,
            "pars3"=>pars3, "traj3"=>traj3, "cost3"=>cost3, "cpm_traj3"=>cpm_traj3, "ftraj3"=>ftraj3,
            "cost"=>cost, "cost1s"=>cost1s, "cost2s"=>cost2s, "bbox"=>bbox, 
            "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))

    catch y
        # Interrupts should not get caught:
        if isa(y, InterruptException); throw(InterruptException()); end

        # Other errors get caught and a warning is issued but then we run again
        append_to_file(report_file, @sprintf("\n\nWhoopsety, unkown error!\n\n"))
        append_to_file(report_file, @sprintf("Error was :\n%s\n\nTrying new random seed.\n\n", y))
    end

    # Change random seed before next iteration so we don't get stuck in one loop
    extra_pars[:seedrand] = extra_pars[:seedrand]+1
    append_to_file(report_file, @sprintf("\n\n\Changing to random seed %d\n\n\n", extra_pars[:seedrand]))
end
    




