# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


README_TOP = """

File to try a non-opto regular search; if that doesn't work, try first using few trials
and a target V[1]-V[4] using new_J() to find seed params; and then try again from there.

"""

include("pro_anti.jl")


# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.


"""
cost, costf, dP, dA, hBP, hBA, proValls, antiValls = 
    new_J(nPro, nAnti; pro_target_diff=0.1, anti_target_diff=0.1, 
    opto_conditions = Array{Any}(0,4), plot_condition = 0,
    verbose=false, pre_string="", seedrand=NaN, 
    rule_and_delay_periods = nothing, target_periods = nothing, post_target_periods = nothing,
    plot_list = [], model_params...)

Runs a proAnti network, with a cost that simply asks for a certain fixed, signed, difference
between the Pro units on Pro trials and another fixed, signed, difference between Anti units.
Note that there is nothing about per cent correct here! E.g., 100% of Pro trials have the same target.
Like JJ(), this function can, if desired, run across multiple opto conditions and multiple period 
durations and returns resulting total cost.

(The motivation is to use a search on this new_J() to seed starting param values for the full JJ().)

If rule_and_delay_periods, target_periods, and post_target_periods are not passed, it tries to get them from 
their singular (not plural) versions in model_params, e.g., model_params[:rule_and_delay_period]. NOTE that this
is not what JJ() does for target_period.

# PARAMETERS:

- nPro, nAnti    The number of Pro and the number of Anti trials to run

- pro_target     The target V[1] - V[4] for Pro trials.

- anti_target    The target V[4] - V[1] for Anti trials

- opto_conditions    Each row is the [start_time, end_time, pro_target, anti_target] for an opto condition. 
                The first two columns follow the rules of `parse_opto_times()`: they can be arbitrary Julia 
                expressions involving the terms trial_start, target_start, target_end, and trial_end.

- plot_condition   If non-zero, should be an integer in the range of the rows of opto_conditions, and
                indicates which condition to plot.  A zero means don't plot any of them.

- seedrand      If sets, calls srand() on the value to initialize the random number generator.

- verbose       If true, prints out diagnostic information to the console.

- pre_string    Relevant only under verbose=true, a string that gets printed out before the rest of the verbose info.

- rule_and_delay_periods    Vector, indicating set of rule_and_delay_period lengths to iterate over, 
                            while testing set of opto_periods, etc. on each one. 
                            Deafult is to do a single one, as picked out from model_params[:rule_and_delay_period]

- target_periods            Vector, indicating set of target_perdiod lengths to iterate over, 
                            while testing set of opto_periods, etc. on each one.
                            ** DEFAULT IS TO USE 0.1**, not to pick it out from model_params

- post_target_periods       Vector, indicating set of post_target_period lengths to iterate over, 
                            while testing set of opto_periods, etc. on each one. 
                            Deafult is to do a single one, as picked out from model_params[:post_target_period]
 
- plot_list                 A list of trial numbers to plot in each condition, e.g.  [1:10;]

- model_params              Any remaining keyword-value params are passed as is on to run_ntrials 


, get_value(costf), get_value(dP), get_value(dA), get_value(hBP), get_value(hBA), 
        get_value(proValls), get_value(antiValls)

# RETURNS:

- cost   The net cost, composed of squared error cost (promoting signed V[1]-V[4] differences close to the desired ones).

- costf  A matrix, with cost for each opto x period_length condition

- dP     squared error cost, mean of (V[1]-V[4] - pro_target)^2 on Pro trials
    
- dA     squared error cost, mean of (V[4]-V[1] - anti_target)^2  on Anti trials
    
- hBP    Pro binarized hits, as computed by binarizing (equivalent to theta1->0)
    
- hBA    Anti binarized hits

- proValls   record of full V as a function of time for Pro trials

- antiValls  record of full V as a function of time for Anti trials


"""

function new_J(nPro, nAnti; pro_target_diff=0.1, anti_target_diff=0.1, 
    opto_conditions = Array{Any}(0,4), plot_condition = 0,
    verbose=false, pre_string="", seedrand=NaN, 
    rule_and_delay_periods = nothing, target_periods = nothing, post_target_periods = nothing,
    plot_list = [], model_params...)

    if FDversion() < 0.6
        error("Sorry, new_J() runs only on Julia 0.6 or higher")
    end
    

    # All the variables that we MIGHT choose to differentiate w.r.t. go into this bag -- further down
    # we'll use get_eltype(varbag) to check for any of them being ForwardDiff.Dual.
    # That is how we'll tell whether new matrices should be regular numbers of ForwardDiff.Dual's.
    # *** if you add a new variable you'll want to differentiate w.r.t., it should be added here too ***
    varbag = (pro_target_diff, anti_target_diff, opto_conditions, model_params)
    # print("get_eltype(varbag)="); print(get_eltype(varbag)); print("\n")
    
    # If the plurals of the periods are not passed in, then use the singular in model_params as the default:
    if rule_and_delay_periods==nothing
        rule_and_delay_periods = model_params[:rule_and_delay_period]
    end
    if target_periods==nothing        
        target_periods = model_params[:target_period]
    end
    if post_target_periods==nothing
         post_target_periods = model_params[:post_target_period]
    end
    
    if size(opto_conditions,1)==0  # if there's no opto that is being asked for
        # Then run with only a single opto_period request, with no opto, and control targets as our targets
        opto_conditions = ["trial_start-1" "trial_start-1" pro_target_diff anti_target_diff]
    end
    
    if size(opto_conditions,2) != 4
        try
            # Make sure its rows are 4 cols
            opto_conditions = reshape(opto_conditions, Int64(round(length(opto_periods)/4)), 4) 
        catch
            error("Something is wrong with opto_periods -- it should have 4 columns")
        end
    end
    
    noptos     = size(opto_conditions,1)  # of opto conditions
    nruns_each = length(rule_and_delay_periods)*length(target_periods)*length(post_target_periods)    # runs per opto condition
    nruns      = nruns_each*noptos  # total conditions

    costf = zeros(get_eltype(varbag), noptos, nruns)
    
    dP  = zeros(noptos, nruns_each);   # Pro  diffs
    dA  = zeros(noptos, nruns_each);   # Anti diffs
    hBP = zeros(noptos, nruns_each);   # Pro binarized hits
    hBA = zeros(noptos, nruns_each);   # Anti binarized hits

    proValls         = [];
    antiValls        = [];
    opto_fraction    = [];
    pro_input        = [];
    anti_input       = [];
    
    for nopto=1:noptos # iterate over each opto inactivation period
        # @printf("size(hBP) is %d, %d\n", size(hBP,1), size(hBP,2))

        # reset random number generator for each opto period, so it cant over fit noise samples
        if ~isnan(seedrand); srand(seedrand); end

        n = 0  # n is a counter over all period duration conditions
        totHitsP = totHitsA = totDiffsP = totDiffsA = 0
        for i in rule_and_delay_periods
            for j in target_periods
                for k = post_target_periods
                    n += 1

                    # include this opto inactivation in the parameters to pass on
                    my_params = make_dict(["rule_and_delay_period","target_period","post_target_period"], [i,j,k])
                    my_params = make_dict(["opto_times"], [reshape(opto_conditions[nopto,1:2], 1, 2)], my_params)
                    my_params = merge(Dict(model_params), my_params)  # my_params takes precedence

                    my_plot_list = []
                    if plot_condition == nopto; my_plot_list = plot_list; else myplot_list=[]; end

                    # print("model params is " ); print(model_params); print("\n")
                    proVs, antiVs, proVall, antiVall, opto_fraction,pro_input,anti_input =
                        run_ntrials(nPro, nAnti; plot_list=my_plot_list, my_params...)

                    if length(proValls)==0
                        proValls = zeros(4, size(proVall,2), size(proVall,3), noptos)
                    end
                    if length(antiValls)==0
                        antiValls = zeros(4, size(antiVall,2), size(antiVall,3), noptos)
                    end
                    proValls[:,:,:,nopto]  = get_value(proVall)  # make sure they're not stored as ForwardDiff Duals
                    antiValls[:,:,:,nopto] = get_value(antiVall)
                    
                    diffsP = proVs[1,:]  - proVs[4,:]
                    diffsA = antiVs[4,:] - antiVs[1,:]
                    
                    # set up storage  -- we do get_value() to make sure to from ForwardDiff.Dual into Float64 if necessary
                    dP[nopto, n] = sqrt(mean((get_value(diffsP) - opto_conditions[nopto,3]).^2));
                    dA[nopto, n] = sqrt(mean((get_value(diffsA) - opto_conditions[nopto,4]).^2));
                    hBP[nopto, n] = get_value(sum(proVs[1,:] .>= proVs[4,:,])/nPro);
                    hBA[nopto, n] = get_value(sum(antiVs[4,:] .>  antiVs[1,:,])/nAnti);                    
                    
                    if nPro>0 && nAnti>0
                        # cost can accept ForwardDiff.Dual, so no get_value() for them
                        costf[nopto, n] = (nPro *mean((diffsP - opto_conditions[nopto,3]).^2) +
                                         nAnti*mean((diffsA - opto_conditions[nopto,4]).^2))/(nPro+nAnti)
                    elseif nPro>0
                        costf[nopto, n] = mean((diffsP - opto_conditions[nopto,3]).^2)
                    else
                        costf[nopto, n] = mean((diffsA - opto_conditions[nopto,4]).^2)
                    end
                end
            end
        end
    
        if verbose
            pcost = mean(costf[nopto,:])   # partial costs
            
            # Notice the get_value() calls below, to transform ForwardDiff Duals into Float64s
            @printf("%s", pre_string)
            @printf("Opto condition # %d\n", nopto)
            @printf("     - %d - cost=%g\n", nopto, get_value(pcost))
            if nPro>0 && nAnti>0
                @printf("     - %d - hBP=%g, dP=%g, hBA=%g, dA=%g\n", nopto, 
                    mean(hBP[nopto,:]), mean(dP[nopto,:]), mean(hBA[nopto,:]), mean(dA[nopto,:]))
            elseif nPro>0
                @printf("     - %d - hBP=%g, dP=%g\n", nopto, mean(hBP[nopto,:]), mean(dP[nopto,:]))
            else
                @printf("     - %d - hBA=%g, dA=%g\n", nopto, mean(hBA[nopto,:]), mean(dA[nopto,:]))
            end        
        end
    end
        
    cost = mean(costf)

    if verbose
        @printf("%s", pre_string)
        @printf("OVERALL\n")
        @printf("     -- cost=%g\n", get_value(cost))
    end
    

    # The scalar cost should be differentiable, the others should be regular Float64s
    return cost, get_value(costf), get_value(dP), get_value(dA), get_value(hBP), get_value(hBA), 
        get_value(proValls), get_value(antiValls)
end                    


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




# DON'T MODIFY THIS FILE -- the source is in file Current Carlos Work.ipynb. Look there for further documentation and examples of running the code.

ftraj2 = []; cost2 = [];
extra_pars[:seedrand] = Int64(round(1000*time()))   # 1510782006169 causes lin.alg error but then looks like it'll succeed
srand(extra_pars[:seedrand])

search_range = extra_pars[:search_range]; 

README = """

Farm C9: Like C8, only do a single full search, first new_J() then full search,
no pre-search. But only use the params found by new_J() if they ended up ok

"""

if !isdir("../NewFarms"); mkdir("../NewFarms"); end
fbasename = "../NewFarms/farm_C9_"
# If we wanted a unique identifier per processor run the following line would help:
# if ~isnull(tryparse(Int64, ARGS[1])); fbasename = fbasename * ARGS[1] * "_"; end

@printf("\n\n\nStarting with random seed %d\n\n\n", extra_pars[:seedrand])

while true
    
    @printf("\n\n--- new run ---\n\n")
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
        bbox = merge(bbox, Dict(k => Array{Float64}(search_conditions[k][3])))
    end
    args = Array{String}(args)
    seed = Array{Float64}(seed)


    maxiter1 = 1000;   # for func1, the regular search
    maxiter2 = 2000;  # for func2, the reduced search
    testruns = 10000; # Number of trials for evaluating the results of the model. 10000 is a good number 


    # Make sure to keep the noise frozen over the search, meaning JJ() needs the seedrand parameter
    func1 =  (;params...) -> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false, 
        merge(merge(mypars, extra_pars), Dict(params))...)[1]
    
    # For func2:
    extra_pars[:opto_conditions] = []    
    extra_pars[:plot_condition] = 0

    # new_J() has to be allowed to return all its outputs, so that
    # stopping_func can use some of them
    func2 =  (;params...) -> new_J(40, 40; pre_string="new_J(): ", 
        verbose=false, merge(merge(mypars, extra_pars), Dict(params))...)

    # For func2:
    # This function will get the output of new_J() at each iteration, and will return "true", stopping
    # the minimization, if hBP and hBA are both above a threshold.
    function stopping_func(;cost=0, func_out=[], ignored_extra_params...)
        costf, dP, dA, hBP, hBA = func_out
        return hBP[1]>=0.6 && hBA[1] >= 0.6
    end


    try
        ntries = 2

        extra_pars[:plot_list] = []

        # First with func2 for new_J()
        pars2, traj2, cost2, cpm_traj2, ftraj2 = bbox_Hessian_keyword_minimization(seed, 
            args, bbox, func2, 
            stopping_function = stopping_func, 
            start_eta = 0.1, tol=1e-12, verbose=true, verbose_every=1, maxiter=maxiter2)

        # If this resulted in a place that look good, use it, otherwise start afresh:
        if stopping_func(;cost=cost2, func_out=ftraj2[3,end]); start_pars=pars2; else start_pars=seed; end

        # Then with func1 for JJ()
        pars3, traj3, cost3, cpm_traj3, ftraj3 = bbox_Hessian_keyword_minimization(pars2, 
            args, bbox, func1, 
            start_eta = 0.1, tol=1e-12, 
            verbose=true, verbose_every=1, maxiter=maxiter1)
            
        cost, cost1s, cost2s, hP, hA, dP, dA, hBP, hBA = JJ(testruns, testruns; verbose=false, 
            make_dict(args, pars3, merge(merge(mypars, extra_pars)))...)
        
        myfilename = next_file(fbasename, 4)
        myfilename = myfilename*".jld"

        @printf("\n\n ****** writing to file %s *******\n\n", myfilename)
        
        # write file
        save(myfilename, Dict("README"=>README, "nPro"=>mypars[:nPro], "nAnti"=>mypars[:nAnti], "ntries"=>ntries, 
            "mypars"=>mypars, "extra_pars"=>extra_pars, "args"=>args, "seed"=>seed, "bbox"=>bbox, 
            "pars2"=>pars2, "traj2"=>traj2, "cost2"=>cost2, "cpm_traj2"=>cpm_traj2, 
            "pars3"=>pars3, "traj3"=>traj3, "cost3"=>cost3, "cpm_traj3"=>cpm_traj3, "ftraj3"=>ftraj3,
            "cost"=>cost, "cost1s"=>cost1s, "cost2s"=>cost2s,
            "hP"=>hP, "hA"=>hA, "dP"=>dP, "dA"=>dA, "hBP"=>hBP, "hBA"=>hBA))

    catch y
        if isa(y, InterruptException); throw(InterruptException()); end
        @printf("\n\nWhoopsety, unkown error!\n\n");
        @printf("Error was :\n"); print(y); @printf("\n\nTrying new random seed.\n\n")
    end

    # Change random seed so we don't get stuck in one loop
    extra_pars[:seedrand] = extra_pars[:seedrand]+1
end


