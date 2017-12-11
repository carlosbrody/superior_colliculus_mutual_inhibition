

include("pro_anti.jl")

mypars = Dict(
:vW     =>                   0, 
:hW     =>                   0,
:dW     =>                   0, 
:sW     =>                   0, 
:sigma  =>                   0,
:constant_excitation      => 0, 
:target_period_excitation => 0, 
:right_light_excitation   => 0, 
:const_pro_bias           => 0, 
:opto_strength            => 0, 
:pro_rule_strength        => 0, 
:anti_rule_strength       => 0, 
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

binarized_delta_threshold = 0.05   # maximum average difference between binarized hits and their targets for a run to go into full optimization

README = """
Optimizing existing C17 farms-- running optimizations, 1600 trials per condition, starting 
from their ending points but with a fresh random seed.
"""

if ~isnull(tryparse(Int64, ARGS[1])); my_run_number = parse(Int64, ARGS[1]); 
else                                  my_run_number = 1; 
end
if ~isnull(tryparse(Int64, ARGS[2])); tot_n_runs    = parse(Int64, ARGS[2]); 
else                                  tot_n_runs = 1; 
end

sourcedir = "../Farms" * chomp(readstring(`hostname`))[end-2:end]
if !isdir(farmdir); mkdir(sourcedir); end
fbasename = sourcedir * "/" * "farm_C17_"

optim_dir = "../Farms" * chomp(readstring(`hostname`))[end-2:end] * "_Optim"
if !isdir(optim_dir); mkdir(optim_dir); end

reports_dir = "../Reports" * chomp(readstring(`hostname`))[end-2:end]
if !isdir(reports_dir); mkdir(reports_dir); end
report_file = reports_dir * "/" * @sprintf("report_out_%d", my_run_number)

f = readdir(source_dir)
nloops = 0

while my_run_number + (nloops*tot_n_runs) <= length(f)
    myfile = f[my_run_number + (nloops*tot_n_runs)]

    append_to_file(report_file, @sprintf("\n\n*** %d %s: probing file %s ***\n\n", my_run_number, 
        Dates.format(now(), "e, dd u yyyy HH:MM:SS"), myfile))

    mypars, extra_pars, hBP, hBA = load(source_dir * "/" * myfile, "mypars", "extra_pars", "hBP", "hBA")
    opto_targets = extra_pars[:opto_targets]
    
    if mean(opto_targets - [hBP hBA]) < binarized_delta_threshold

        
        append_to_file(report_file, @sprintf("\n\n*** %d %s: grabbing file %s ***\n\n", my_run_number, 
                                             Dates.format(now(), "e, dd u yyyy HH:MM:SS"), myfile))

        mypars, extra_pars, args, seed = load(source_dir * "/" * myfile, 
                                              "mypars", "extra_pars", "args", "pars3")
        mypars[:nPro]  = 1600
        mypars[:nAnti] = 1600

        append_to_file(report_file, @sprintf("Running with %d Pro trials and %d Anti trials.\n\n", 
                                             mypars[:nPro], mypars[:nAnti]))

        extra_pars[:seedrand] = Int64(round(time()*1000))

        maxiter1 = 1000;   # for func1, the regular search
        testruns = 10000;  # Number of trials for evaluating the results of the model. 10000 is a good number 

        func1 =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, 
                                    merge(merge(mypars, extra_pars), Dict(params))...)[1]
        try
            pars3, traj3, cost3, cpm_traj3, ftraj3 = 
                bbox_Hessian_keyword_minimization(seed, 
                                                  args, bbox, func1, 
                                                  start_eta = 0.01, tol=1e-12, 
                                                  verbose_file = report_file,
                                                  verbose=true, verbose_every=1, maxiter=maxiter1)

            # evaluate the result with many trials, for accuracy
            cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(testruns, testruns; verbose=false, 
                                                                make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)

            # Write out the results
            myfilename = optim_dir * "/" * myfile

            append_to_file(report_file, @sprintf("\n\n ****** writing to file %s *******\n\n", myfilename))

            # write file
            save(myfilename, Dict("README"=>README, "nPro"=>mypars[:nPro], "nAnti"=>mypars[:nAnti], 
                                  "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, 
                                  "pars3"=>pars3, "traj3"=>traj3, "cost3"=>cost3, "cpm_traj3"=>cpm_traj3, "ftraj3"=>ftraj3,
                                  "cost"=>cost, "cost1s"=>cost1s, "cost2s"=>cost2s, "bbox"=>bbox,
                                  "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))

        catch y
            # Interrupts should not get caught:
            if isa(y, InterruptException); throw(InterruptException()); end

            # Other errors get caught and a warning is issued but then we run again
            append_to_file(report_file, @sprintf("\n\nWhoopsety, unkown error!\n\n"));
            append_to_file(report_file, @sprintf("Error was:\n%s\n\nTrying new random seed.\n\n", y)); 
        end
    end
    nloops += 1
end


