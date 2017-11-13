# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


include("pro_anti.jl")

#
#  We're assuming that results of reduced_farm.jl were already generated.
#


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


"""
function show_run(run_name, ntrials=100; further_pars...)

Given a string representing the filename of a run generated in reduced_farm.jl
(for example, "farm_C3_0027.jld"), tries to load that file from within "FarmFields/"
runs it to plot 30 Pro trials in figure 1 and 30 Anti trials in figure 2,
and indicates in those figures the corresponding number of binarized hits.

Any further parameters override the contents of the farm and are passed on to 
run_ntrials(), so for example you can set rule_and_delay_period=0.2
"""
function show_run(run_name, ntrials=100; further_pars...)

    pygui(true)
    figure(1); clf();
    figure(2); clf();

    d = load("FarmFields/" * run_name)

    extra_pars = d["extra_pars"]
    mypars     = d["mypars"]
    args       = d["args"]
    pars       = d["pars"]

    delete!(extra_pars, :plot_list)   # we want our own plot_list
    delete!(extra_pars, :seedrand)    # don't use the farm's seedrand. (Could change this if wanted)

    # Add or replace alternative period values if needed
    # if rule_and_delay_period != nothing; extra_pars[:rule_and_delay_period] = rule_and_delay_period; end;
    # if target_period         != nothing; extra_pars[:target_period]         = target_period;         end;
    # if post_target_period    != nothing; extra_pars[:post_target_period]    = post_target_period;    end;

    # Run it
    proVs, antiVs = run_ntrials(ntrials, ntrials; plot_list=[1:30;], 
        merge(make_dict(args, pars, merge(mypars, extra_pars)), Dict(further_pars))...)[1:2]

    # Add titles with the binarized hit rates
    hBP = 100*length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
    hBA = 100*length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)

    figure(1); subplot(3,1,1); title(@sprintf("%s: binarized Pro  hits = %.1f %%\n", run_name, hBP))
    figure(2); subplot(3,1,1); title(@sprintf("%s: binarized Anti hits = %.1f %%\n", run_name, hBA))
    
    return hBP, hBA
end



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# IDENTIFY GOOD RUNS FROM C3 FARM

I = results = files = 0

recompute_me = false; if recompute_me
    results = zeros(0,4)
    files = []
    for f in filter(x -> startswith(x, "farm_C3_"), readdir("FarmFields/"))
        hBP, hBA, dP, dA= load("FarmFields/" * f, "hBP", "hBA", "dP", "dA")
        results = [results; [hBP hBA dP dA]]; 
        files   = [files ; f]
    end
    
    I = find((results[:,1] .> 0.9) .& (results[:,2] .> 0.65) .& (results[:,3].>0.85) .& (results[:,4].>0.85))

    @save "Temp/C3_results.jld" results files I
else
    I, results, files = load("Temp/C3_results.jld", "I", "results", "files");
end



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


"""
    go_over_good_runs(files, I)

Given a list of files, and a list of indices into those files, iterates over the indicated
filenames, running show_run() on them, and then asking for user input.
"""
function go_over_good_runs(files, I)
    for k=1:length(I)
        show_run(files[I[k]])

        @printf("Any key to continue, q to quit: ")
        ans = chomp(readline())
        if ans=="q" || ans=="Q"
            break
        end
    end
end



# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


# show_run("farm_C3_0027.jld", 400; opto_times=["trial_start", "trial_end"], opto_strength=0.75)
# show_run("farm_C3_0027.jld", 400; opto_times=["target_start-0.5", "target_start"], opto_strength=0.4, 
#    const_pro_bias=0.02476, anti_rule_strength=0.05)
show_run("farm_C3_0027.jld", 400; opto_times=["target_start+0.15", "target_end"], opto_strength=0.4, 
    const_pro_bias=0.02476, anti_rule_strength=0.054)


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


#  TARGET PERIOD OPTO

run_name = "farm_C3_0027.jld"
recompute_me = false; if recompute_me
    ostrengths = 0.5:0.1:1
    phits      = zeros(1, length(ostrengths))
    ahits      = zeros(1, length(ostrengths))

    for i=1:length(ostrengths)
        phits[i], ahits[i] = show_run(run_name, 4000; 
        opto_times=["target_start+0.016", "target_end+0.016"], opto_strength=ostrengths[i])
        @printf("ostrength=%g, phits=%.1f, ahits=%.1f\n", ostrengths[i], phits[i], ahits[i])
    end

    @save("Temp/opto_target.jld", ostrengths, ahits, phits)
else
    ostrengths, ahits, phits = load("Temp/opto_target.jld", "ostrengths", "ahits", "phits")
end

figure(3);clf();
plot(ostrengths, phits[:], ".-", ostrengths, ahits[:], ".-")
legend(["pro", "anti"])
xlabel("opto strength")
ylabel("binarized hits")
title(run_name * " : TARGET PERIOD opto")
grid("on")


#  DELAY PERIOD OPTO  ---------

recompute_me = false; if recompute_me
    ostrengths = 0.5:0.1:1
    phits      = zeros(1, length(ostrengths))
    ahits      = zeros(1, length(ostrengths))

    for i=1:length(ostrengths)
        phits[i], ahits[i] = show_run(run_name, 4000; 
        opto_times=["target_start-0.5", "target_start"], opto_strength=ostrengths[i])
        @printf("ostrength=%g, phits=%.1f, ahits=%.1f\n", ostrengths[i], phits[i], ahits[i])
    end
    @save("Temp/opto_delay.jld", ostrengths, ahits, phits)
else
    ostrengths, ahits, phits = load("Temp/opto_delay.jld", "ostrengths", "ahits", "phits")
end

figure(4);clf();
plot(ostrengths, phits[:], ".-", ostrengths, ahits[:], ".-")
legend(["pro", "anti"])
xlabel("opto strength")
ylabel("binarized hits")
title(run_name * " : last 500 ms of DELAY PERIOD opto")
grid("on")


