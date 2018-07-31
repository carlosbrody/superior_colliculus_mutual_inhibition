include("setup_C31_control.jl")


README = """

Farm C31: 3 node network, similar to Farm C30. CONTROL ONLY

"""

if length(ARGS)>0  &&  ~isnull(tryparse(Int64, ARGS[1])); 
    my_run_number = parse(Int64, ARGS[1]); 
else                                                      
    my_run_number = 1; # I am process my_run_number
end

if length(ARGS)>1  &&  ~isnull(tryparse(Int64, ARGS[2])); 
    tot_n_runs    = parse(Int64, ARGS[2]); 
else
    tot_n_runs = 1;   # I'm being run as part of tot_n_run processes
end


farmdir = "../Farms_C31_control"
if !isdir(farmdir); mkdir(farmdir); end
fbasename = farmdir * "/" * "farm_C31_control_" * chomp(readstring(`hostname`)) * "_"

reports_dir = "../Reports" 
if !isdir(reports_dir); mkdir(reports_dir); end
report_file = reports_dir * "/" * @sprintf("report_%s_%d", chomp(readstring(`hostname`)), my_run_number)
if isfile(report_file); rm(report_file); end

extra_pars[:seedrand]                  = Int64(round(time()*parse(Int64,ARGS[1])))
extra_pars[:maxiter]                   = 1000
extra_pars[:testruns]                  = 10000
extra_pars[:few_trials]                = 50
extra_pars[:many_trials]               = 1600
extra_pars[:binarized_delta_threshold] = 0.1    # average frac correct must be within this of target
extra_pars[:anti_perf_delta]           = 0.05   # delay anti must be at least this worse off than control or choice anti
extra_pars[:pro_better_than_anti]      = true   # if true, in each condition pro hits must be > anti hits

append_to_file(report_file, @sprintf("\n\n\nStarting with random seed %d\n\n\n", extra_pars[:seedrand]))

just_testing = false; if just_testing   # 
    extra_pars[:maxiter]                   = 4
    extra_pars[:testruns]                  = 100
    extra_pars[:few_trials]                = 10
    extra_pars[:many_trials]               = 100
    extra_pars[:binarized_delta_threshold] = 0.5     # average frac correct must be within this of target
    extra_pars[:anti_perf_delta]           = 0       # delay anti must be at least this worse off than control or choice anti
    extra_pars[:pro_better_than_anti]      = false   # if true, in each condition pro hits must be > anti hits
end


#while true
    
    append_to_file(report_file, @sprintf("\n\n--- new run -- %s ---\n\n", Dates.format(now(), "e, dd u yyyy HH:MM:SS")))

    args = []; seed = []; bbox = Dict()
    for k in keys(search_conditions)
        search_box = search_conditions[k][2]
        args = [args; String(k)]
        # --- search within the full indicated search box
        seed = [seed ; rand()*diff(search_box) + search_box[1]]
        bbox = merge(bbox, Dict(k => Array{Float64}(search_conditions[k][3])))
    end
    args = Array{String}(args)
    seed = Array{Float64}(seed)

    extra_pars[:nPro]     = extra_pars[:few_trials]
    extra_pars[:nAnti]    = extra_pars[:few_trials]

    # Make sure to keep the noise frozen over the search, meaning JJ_3node() needs the seedrand parameter
    func_quiet =  (;params...) -> JJ_3node(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false, 
        merge(merge(mypars, extra_pars), Dict(params))...)[1]
    
    # try
        parsA, trajA, costA, cpm_trajA, ftrajA = bbox_Hessian_keyword_minimization(seed, 
            args, bbox, func_quiet, 
            start_eta = 0.01, tol=1e-12, verbose_file=report_file, verbose_timestamp=true,
            verbose=true, verbose_every=1, maxiter=extra_pars[:maxiter])
            
        # evaluate the result with many trials, for accuracy
        costA, cost1sA, cost2sA, hPA, hAA, dPA, dAA, hBPA, hBAA = JJ_3node(extra_pars[:testruns], extra_pars[:testruns]; 
                                                   verbose=false, 
                                                   make_dict(args, parsA, merge(merge(mypars, extra_pars)))...)
        
        append_to_file(report_file, @sprintf("mean diff = %g\n", mean(abs.(extra_pars[:opto_targets] - [hBPA hBAA]))))
        append_to_file(report_file, @sprintf("all(hBPA .> hBAA) = %s\n", all(hBPA .> hBAA)))
        print_vector_g(report_file, hBPA); print_vector_g(report_file, hBAA)
        append_to_file(report_file, @sprintf("\nhBAA[2] <= hBAA[1] - extra_pars[:anti_perf_delta] = %s\n", 
                                             hBAA[2] <= hBAA[1] - extra_pars[:anti_perf_delta]))
        append_to_file(report_file, @sprintf("hBAA[2] <= hBAA[3] - extra_pars[:anti_perf_delta] = %s\n", 
                                             hBAA[2] <= hBAA[3] - extra_pars[:anti_perf_delta]))


        # always completes further training

            append_to_file(report_file, @sprintf("\n\n**** training further **** %s ---\n\n", 
                                                 Dates.format(now(), "e, dd u yyyy HH:MM:SS")))

            extra_pars[:nPro]     = extra_pars[:many_trials]
            extra_pars[:nAnti]    = extra_pars[:many_trials]

            func_chatty =  (;params...) -> JJ_3node(extra_pars[:nPro], extra_pars[:nAnti]; verbose=true,
                                              verbose_file = report_file,
                                              merge(merge(mypars, extra_pars), Dict(params))...)[1]
            try

                pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(parsA, 
                                                                 args, bbox, func_chatty,  verbose_timestamp=true,
                                                                 start_eta = 0.01, tol=1e-12, verbose_file=report_file,
                                                                 verbose=true, verbose_every=1, maxiter=extra_pars[:maxiter])
            
                # evaluate the result with many trials, for accuracy
                cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ_3node(extra_pars[:testruns], extra_pars[:testruns]; 
                                                               verbose=false, 
                                                               make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)

                myfilename = next_file(fbasename, 4)
                myfilename = myfilename*".jld"

                append_to_file(report_file, @sprintf("\n\n ****** writing to file %s *******\n\n", myfilename))
        
                # write file
                save(myfilename, Dict("README"=>README, "nPro"=>extra_pars[:nPro], "nAnti"=>extra_pars[:nAnti], 
                                      "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, "bbox"=>bbox, 
                                      "search_conditions"=>search_conditions,
                                      "parsA"=>parsA, "trajA"=>trajA, "costA"=>costA, "cpm_trajA"=>cpm_trajA, 
                                      "ftrajA"=>ftrajA,
                                      "costA"=>costA, "cost1sA"=>cost1sA, "cost2sA"=>cost2sA, 
                                      "hPA"=>hPA, "hAA"=>hAA, "dPA"=>dPA, "dAA"=>dAA, "hBPA"=>hBPA, "hBAA"=>hBAA,            
                                      "pars3"=>pars3, "traj3"=>traj3, "cost3"=>cost3, "cpm_traj3"=>cpm_traj3, 
                                      "ftraj3"=>ftraj3,
                                      "cost"=>cost, "cost1s"=>cost1s, "cost2s"=>cost2s, "bbox"=>bbox, 
                                      "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))

            catch y1
                # Interrupts should not get caught:
                if isa(y1, InterruptException); throw(InterruptException()); end

                # Other errors get caught and a warning is issued but then we run again
                append_to_file(report_file, @sprintf("\n\nWhoopsety, unkown error during LONG search!\n\n"))
                append_to_file(report_file, @sprintf("Error was :\n%s\n\nTrying new random seed.\n\n", y1))
            end
#    catch y
        # Interrupts should not get caught:
#        if isa(y, InterruptException); throw(InterruptException()); end

        # Other errors get caught and a warning is issued but then we run again
#        append_to_file(report_file, @sprintf("\n\nWhoopsety, unkown error during QUICK search!\n\n"))
#        append_to_file(report_file, @sprintf("Error was :\n%s\n\nTrying new random seed.\n\n", y))
#    end
    
    extra_pars[:seedrand] = extra_pars[:seedrand]+1
    append_to_file(report_file, @sprintf("\n\n\Changing to random seed %d\n\n\n", extra_pars[:seedrand]))
#end
    




