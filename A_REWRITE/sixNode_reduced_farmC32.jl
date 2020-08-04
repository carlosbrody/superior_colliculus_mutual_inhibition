include("startup.jl")
include("sixNodeSetup_C32.jl")


README = """

Six Node C32: Just like C32, but with two Anti nodes per side of the brain

"""

s = ArgParseSettings()
@add_arg_table! s begin
    "--report_file"
        help = "the string representing where the output will be logged"
        arg_type = String
    "--bspockID", "-b"
        help = "the job ID that the spock respawner is keeping track of for iterations of this process"
        arg_type = Int64
    "--onSpock", "-s"
        help = "flag used to indicate we're running on spock"
        action = :store_true
    "--respawn", "-r"
        help = "we should respawn, using the logfile for initial conditions. "*
            "Otherwise start up from scratch."
        action = :store_true
    "--hostname"
        help = "string used within farm and report filenames"
        arg_type = String
        default = String(chomp(read(`hostname`, String)))
    "my_run_number"
        help = "if not on spock, will be first positional argument"
        arg_type = Int64
        default = 1
    "tot_n_runs"
        help = "if not on spock, will be second positional argument"
        arg_type = Int64
        default = 1
end
parsed_args = parse_args(myargs, s)
for k in keys(parsed_args)
    eval(Meta.parse("$k = parsed_args[\"$k\"]"))
end
if onSpock
    hostname = "spock"
end


farmdir = "../../Farms_6Node"
if !isdir(farmdir); mkdir(farmdir); end
fbasename = "$farmdir/farm_6Node_$(hostname)_"

reports_dir = "../../Reports"
if !isdir(reports_dir); mkdir(reports_dir); end
if report_file == nothing
    if onSpock
        report_file = "$reports_dir/r6_$(hostname)_$bspockID"
    else
        report_file = "$reports_dir/r6_$(hostname)_$my_run_number"
    end
end

extra_pars[:maxiter]                   = 1000
extra_pars[:testruns]                  = 10000
extra_pars[:few_trials]                = 50
extra_pars[:many_trials]               = 1600
extra_pars[:binarized_delta_threshold] = 0.1    # average frac correct must be within this of target
extra_pars[:anti_perf_delta]           = 0.05   # delay anti must be at least this worse off than control or choice anti
extra_pars[:pro_better_than_anti]      = true   # if true, in each condition pro hits must be > anti hits

just_testing = false; if just_testing   #
    extra_pars[:maxiter]                   = 25
    extra_pars[:testruns]                  = 3
    extra_pars[:few_trials]                = 10
    extra_pars[:many_trials]               = 20
    extra_pars[:binarized_delta_threshold] = 1     # average frac correct must be within this of target
    extra_pars[:anti_perf_delta]           = 0     # delay anti must be at least this worse off than control or choice anti
    extra_pars[:pro_better_than_anti]      = false   # if true, in each condition pro hits must be > anti hits
end

second_stage_training = false
if !respawn
    if isfile(report_file); rm(report_file); end

    extra_pars[:seedrand]                  = Int64(my_run_number*round(time()*1000))
    append_to_file(report_file, "\n\n\nStarting with random seed $(extra_pars[:seedrand])\n\n\n")
    start_eta = 0.01
else
    fp = open(report_file, "r")
    s = read(fp, String)
    close(fp)
    extra_pars[:seedrand] = parseLast(s, "random seed ", Int64)
    start_eta             = parseLast(s, " eta=", Float64)
    seed                  = parseLast(s, " ps=", Vector{Float64})

    tf = findlast("training further", s)
    if tf != nothing
        if tf[1] > findlast("random seed", s)[1]
            second_stage_training = true
        end
    end
end

# --- define the args string vector and bbox for seacrh space bounding box
args = []; bbox = Dict()
for k in keys(search_conditions)
    search_box = search_conditions[k][2]
    global args = [args; String(k)]
    # --- search within the full indicated search box
    global bbox = merge(bbox, Dict(k => Array{Float64}(search_conditions[k][3])))
end
args = Array{String}(args)

while true
    global second_stage_training, respawn, report_file, bspockID, onSpock
    global hostname, seed

    if !second_stage_training
        if !respawn
            append_to_file(report_file, "\n\n--- new run -- "*
                "$(Dates.format(now(), "e, dd u yyyy HH:MM:SS")) ---\n\n")

            # Choose a random starting point:
            seed = [];
            for k in keys(search_conditions)
                # --- search within the full indicated search box
                search_box = search_conditions[k][2]
                seed = [seed ; rand()*diff(search_box) .+ search_box[1]]
            end
            seed = Array{Float64}(seed)
        end

        extra_pars[:nPro]     = extra_pars[:few_trials]
        extra_pars[:nAnti]    = extra_pars[:few_trials]

        # Make sure to keep the noise frozen over the search, meaning JJ() needs the seedrand parameter
        func_quiet =  (;params...) -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=false,
            merge(merge(mypars, extra_pars), Dict(params))...)[1]

        # try
        parsA, trajA, costA, cpm_trajA, ftrajA = bbox_Hessian_keyword_minimization(seed,
            args, bbox, func_quiet,
            start_eta = start_eta, tol=1e-12, verbose_file=report_file, verbose_timestamp=true,
            verbose=true, verbose_every=10, maxiter=extra_pars[:maxiter])

        # evaluate the result with many trials, for accuracy
        costA, cost1sA, cost2sA, hPA, hAA, dPA, dAA, hBPA, hBAA =
            JJ(extra_pars[:testruns], extra_pars[:testruns]; verbose=false,
                make_dict(args, parsA, merge(merge(mypars, extra_pars)))...)

        ## fix all variables, is mean the best, or should I expand the test so each duration works?
        if size(hBPA,2) > 1
            println(report_file, hBPA);
            println(report_file, hBAA);
            hBPA = mean(hBPA, dims=2);
            hBAA = mean(hBAA, dims=2);
        end

        append_to_file(report_file,
            "mean diff = $(mean(abs.(extra_pars[:opto_targets] - [hBPA hBAA])))\n")
        append_to_file(report_file, "all(hBPA .> hBAA) = $(all(hBPA .> hBAA))\n")
        println(report_file, hBPA); println(report_file, hBAA)
        append_to_file(report_file,
            "\nhBAA[2] <= hBAA[1] - extra_pars[:anti_perf_delta] "*
            "= $(hBAA[2] <= hBAA[1] - extra_pars[:anti_perf_delta])\n")
        append_to_file(report_file,
            "hBAA[2] <= hBAA[3] - extra_pars[:anti_perf_delta] = "*
            "$(hBAA[2] <= hBAA[3] - extra_pars[:anti_perf_delta])\n")

        if mean(abs.(extra_pars[:opto_targets] - [hBPA hBAA])) < extra_pars[:binarized_delta_threshold]  &&
            (!extra_pars[:pro_better_than_anti] || all(hBPA .> hBAA))  &&
            (hBAA[2] <= hBAA[1] - extra_pars[:anti_perf_delta]) &&
            (hBAA[2] <= hBAA[3] - extra_pars[:anti_perf_delta]) &&
            hBAA[1] > 0.5  &&
            hBAA[3] > 0.5

            append_to_file(report_file, "\n\n**** training further **** " *
                "$(Dates.format(now(), "e, dd u yyyy HH:MM:SS")) ---\n\n")

            second_stage_training = true
            seed = parsA
        end
    end

    if second_stage_training
        extra_pars[:nPro]     = extra_pars[:many_trials]
        extra_pars[:nAnti]    = extra_pars[:many_trials]

        func_chatty =  (;params...) -> JJ(extra_pars[:nPro], extra_pars[:nAnti]; verbose=true,
                                              verbose_file = report_file,
                                              merge(merge(mypars, extra_pars), Dict(params))...)[1]
            try
                pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(seed,
                     args, bbox, func_chatty,  verbose_timestamp=true,
                     start_eta = start_eta, tol=1e-12, verbose_file=report_file,
                     verbose=true, verbose_every=10, maxiter=extra_pars[:maxiter])

                # evaluate the result with many trials, for accuracy
                cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA =
                    JJ(extra_pars[:testruns], extra_pars[:testruns]; verbose=false,
                        make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)

                myfilename = next_file(fbasename, 4)
                myfilename = myfilename*".jld"

                append_to_file(report_file, "\n\n ****** writing to file $myfilename *******\n\n")

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
                append_to_file(report_file, "\n\nWhoopsety, unkown error during LONG search!\n\n")
                append_to_file(report_file, "Error was :\n$y1\n\nTrying new random seed.\n\n")
            end
        end
#    catch y
        # Interrupts should not get caught:
#        if isa(y, InterruptException); throw(InterruptException()); end

        # Other errors get caught and a warning is issued but then we run again
#        append_to_file(report_file, @sprintf("\n\nWhoopsety, unkown error during QUICK search!\n\n"))
#        append_to_file(report_file, @sprintf("Error was :\n%s\n\nTrying new random seed.\n\n", y))
#    end

    extra_pars[:seedrand] = extra_pars[:seedrand]+1
    append_to_file(report_file, "\n\nChanging to random seed $(extra_pars[:seedrand])\n\n\n")

    second_stage_training = false
    respawn = false
end
