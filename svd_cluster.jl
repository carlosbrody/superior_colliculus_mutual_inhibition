# Set of functions for doing SVD clustering on neural trajectories
#
# To perform analysis on a new farm;
# > results = load_farm_params(farm_id)
# > response = build_response_matrix(farm_id)
# > SVD_analysis(farm_id)

# I have no idea why I need to include these, but load() won't work unless I do!
include("pro_anti.jl")
using HDF5

"""
    load_farm_params(farm_id; farmdir="MiniFarms", verbose=true, verbose_every=50)

Load final parameters for each farm run in <MiniFarms/farm_id>. Also loads the training and test cost for each farm. Saves a file <farm_id>_results.jld 

# OPTIONAL PARAMETERS:

- farmdir   Direction where farm runs are located

- verbose   If true, displays some information while running

- verbose_every How often to display information

# RETURNS

A dict with fields "dirs", "files", "tcost", "cost", and "params"

"""
function load_farm_params(farm_id; farmdir="MiniFarms", verbose=true, verbose_every=50);
    if typeof(farmdir)==String; farmdir=[farmdir]; end
    results = Dict(); dirs=[]; files =[]; qs=[]; tcosts=[]; costs=[]; pars=[]; n=0;
    for dd in farmdir
        for f in filter(x -> startswith(x, "farm_" * farm_id * "_"), readdir(dd * "/"))
            n +=1
            myfile = dd * "/" * f; 
            args, params, traj3, cost = load(myfile, "args", "pars3", "traj3", "cost")
            files = [files ; myfile]; dirs = [dirs ; dd]
            tcosts = [tcosts; traj3[2,end]]; costs = [costs; cost]
            if length(pars)==0; pars=params'; else pars = [pars ; params']; end
            if verbose && rem(n, verbose_every)==0
                @printf("%s %g\n", myfile, tcosts[end])
            end

        end
    end

    results["dirs"] = dirs
    results["files"]  = files
    results["tcost"]  = tcosts
    results["cost"]   = costs
    results["params"] = pars

    myfilename = farmdir[1]*farm_id*"_results.jld";
    save(myfilename, results)
    return results
end


"""
    run_farm(filename; testruns=200, overrideDict=Dict(),all_conditions=false)

Runs the farm at <filename>, <testruns> times. Then computes the average trajectory for each condition: hits/errors X pro/anti. 

# OPTIONAL PARAMETERS:

- testruns Number of runs to compute

- overrideDict  If you want to override some model parameters

- all_conditions If true, returns average trajectory for each opto condition as well 
# RETURNS

a row vector that contains the average trajectory in each trial type and each opto condition

"""
function run_farm(filename; testruns=200, overrideDict=Dict(),all_conditions=false)
    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")
    numConditions = 1;
    if all_conditions
        numConditions = size(extra_pars[:opto_periods],1);
    end
    avgHP = [];
    avgMP = [];
    avgHA = [];
    avgMA = [];
    farm_response = [];
    for period = 1:numConditions
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        ))

        proVs, antiVs, proF, antiF = run_ntrials(testruns, testruns; merge(make_dict(args, pars3, these_pars), overrideDict)...);

        hitBP = find(proVs[1,:]  .> proVs[4,:])
        avghp = mean(proF[:,:,hitBP],3);

        missBP= find(proVs[1,:]  .< proVs[4,:])
        avgmp = mean(proF[:,:,missBP],3);

        hitBA = find(antiVs[4,:] .> antiVs[1,:])
        avgha = mean(antiF[:,:,hitBA],3);

        missBA= find(antiVs[4,:] .< antiVs[1,:])              
        avgma = mean(antiF[:,:,missBA],3);
    
        avgHP = [avghp[1,:,1]; avghp[2,:,1]; avghp[3,:,1]; avghp[4,:,1]];
        avgMP = [avgmp[1,:,1]; avgmp[2,:,1]; avgmp[3,:,1]; avgmp[4,:,1]];
        avgHA = [avgha[1,:,1]; avgha[2,:,1]; avgha[3,:,1]; avgha[4,:,1]];
        avgMA = [avgma[1,:,1]; avgma[2,:,1]; avgma[3,:,1]; avgma[4,:,1]];
        if period == 1
            farm_response = [avgHP' avgMP' avgHA' avgMA'];
        else
            farm_response = [farm_response avgHP' avgMP' avgHA' avgMA'];
        end
    end
#    return avgHP, avgMP, avgHA, avgMA 
    return farm_response, numConditions
end

"""
    build_response_matrix(farm_id)
    
For each farm run, computes the average trajectory for hits/errors X pro/anti. Arranges these average trajectories in a #runs X length(hits/errors X pro/anti)*#timepoints matrix. Saves a file <farm_id>_SVD_response_matrix.jld 

# OPTIONAL PARAMETERS:

- farmdir   Direction where farm runs are located
    
- all_conditions If true, compute response matrix for all opto conditions

# RETURNS

Response matrix

"""
function build_response_matrix(farm_id; farmdir="MiniFarms", all_conditions = false)

    results = load(farmdir*farm_id*"_results.jld");
    response_matrix = [];
    numConditions = 1;
    for i = 1:length(results["cost"])
        filename = results["files"][i];
        farm_response, numC = run_farm(filename;all_conditions=all_conditions);
        numConditions = numC;

        if isempty(response_matrix)
            response_matrix = farm_response;
        else
            response_matrix = [response_matrix; farm_response];
        end
        @printf("%g %s %g\n",i, filename, results["tcost"][i])
    end

    if all_conditions & (numConditions > 1)
        myfilename = farmdir*farm_id*"_SVD_response_matrix"*string(numConditions)*".jld";
    else
        myfilename = farmdir*farm_id*"_SVD_response_matrix.jld";
    end
    save(myfilename, Dict("response"=>response_matrix, "results"=>results, "numConditions"=>numConditions ))

    return response_matrix
end


"""
    SVD_analysis(farm_id; opto_conditions = 1)    

Plots a series of analyses based on the SVD of the average neural response

# PARAMETERS:

- farm_id   Which farm to analyze   

 # OPTIONAL PARAMETERS:
   
- opto_conditions Number of opto_conditions in SVD_response_matrix

- compute_good_only compute SVD only on good farms

- threshold for determining good farms

# RETURNS

None

"""
function SVD_analysis(farm_id; farmdir="MiniFarms", opto_conditions = 1, compute_good_only=false, threshold=-0.0002)

    # Load responses from all models
    if opto_conditions > 1
    response, results = load(farmdir*farm_id*"_SVD_response_matrix"*string(opto_conditions)*".jld", "response","results")
    else
    response, results = load(farmdir*farm_id*"_SVD_response_matrix.jld", "response","results")
    end

    # need to filter out NaN rows 
    # Some farms in some conditions have no errors, so we have NaNs
    nanrows = any(isnan(response),2);
    
    # Filter out farms with bad cost
    if compute_good_only
        tcost = copy(results["tcost"]);
        badcost = tcost .>= threshold;
        nanrows = badcost | nanrows;
    end
    
    r_all = response[!vec(nanrows),:];

    # figure out size of response
    nt = size(r_all,2); # length of each row
    ntOptoCondition = Int(nt/opto_conditions); # length of row for each opto_condition
    ntTrialType = Int(ntOptoCondition/4); # length of row for each trial type hit/miss X pro/anti for each opto condition
    ntNode  = Int(ntTrialType/4); # length of row for each node in each trial type.

    # example of data format
    # Sort out responses in each trial type in control trials
    # Each is 4 nodes x 76 timesteps = 304 
    # node 1, 2, 3, 4 is the same ordering as run_nTrials
    r_hp = r_all[:,1:ntTrialType];
    r_mp = r_all[:,ntTrialType*2+1:2*ntTrialType];
    r_ha = r_all[:,ntTrialType*2+1:3*ntTrialType];
    r_ma = r_all[:,ntTrialType*3+1:end];

    # F[:U], F[:S], F[:Vt]
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);

    # do PCA for the hell of it, compute variance explained
    C = cov(r_all);
    vals, vecs = eig(C);
    vtotal = sum(vals)
    vc  = cumsum(vals);
    varexp = vc./vtotal;
    varexp = flipdim(varexp,1);
   
    # compute variance explained from SVD
    S = copy(F[:S]);
    S = S.^2;
    stotal = sum(S);
    S = flipdim(S,1);
    sc = cumsum(S);
    svarexp = sc./stotal;
    svarexp = flipdim(svarexp,1);

    # Cumulative variance explained. Most in < 20 dimensions
    figure()
    plot(1-varexp)
    plot(1-svarexp)
    title("Cumulative Variance explained PCA and SVD")
    xlabel("dim")
    ylabel("cumulative var explained")

    # PCA and SVD have the same variance explained
    figure()
    plot(vals)
    plot(S/(size(r_all,2)))
    title("Variance explained in each dimension")
    xlabel("dim")
    ylabel("var explained")

    # make scatter plot against SVD loadings
    figure()
    u = copy(F[:U])
    scatter(u[:,1],u[:,2])
    title("SVD U columns 1 and 2")
    ylabel("SVD Dim 2")
    xlabel("SVD Dim 1")   
    
    # Same scatter plot with dot size proportional to cost
    tcost = copy(results["tcost"])
    tcost = tcost[!vec(nanrows)];
    tcost = convert(Array{Float64,1}, tcost);
    tcost = -tcost;
    tcost = tcost - minimum(tcost);
    tcost = tcost./maximum(tcost);
    figure()
    scatter(u[:,1],u[:,2],s=tcost)
    title("SVD U columns 1 and 2")
    ylabel("SVD Dim 2")
    xlabel("SVD Dim 1")   
    
    # look at histograms of loadings into each dimension. Dim 1 strongly bimodal
    figure()
    h = plt[:hist](u[:,1],40)
    ylabel("SVD Dim 1")
    xlabel("count")
  
    # Dim 2 unimodal
    figure()
    h = plt[:hist](u[:,2],40)
    ylabel("SVD Dim 2")
    xlabel("count")

    # Dim 3 slightly bimodal  
    figure()
    h = plt[:hist](u[:,3],40)
    ylabel("SVD Dim 3")
    xlabel("count")
    # all higher dimensions strongly unimodal
    
    # Because dim 1 and 3 are bimodal, lets scatter with respect to them
    figure()
    scatter(u[:,1],u[:,3])
    title("SVD U columns 1 and 3")
    # Looks like maybe three clusters?
    ylabel("SVD Dim 3")
    xlabel("SVD Dim 1")   

    # We don't see clustering of parameters, so lets look for what parameters correlate with svd dim 1 and 3
    p = copy(results["params"]);
    p = p[!vec(nanrows),:];
    meanp = mean(p,1);
    stdp  = std(p,1);

    # gotta z-score parameters
    for i=1:size(p,2)
        p[:,i] -= mean(p[:,i]);
        p[:,i] /= std(p[:,i]);
    end
    unorm = copy(u[:,1])
    unorm -= mean(unorm)
    unorm /= std(unorm)
    x = [unorm p];
    # compute covariance between each parameter and svd dim 1 
    C = cov(x);
    # Grab parameter labels
    args = load(results["files"][1],"args");
    # Make it a unit vector
    cn = C[1,2:end]/norm(C[1,2:end]);
    # Here is the interesting part!    
    c_labels = [cn args]

    # Lets do the same analysis with svd-dim3
    # u3 is orthogonal to u1, so we dont expect to get the same thing
    unorm3 = copy(u[:,3])
    unorm3 -= mean(unorm3)
    unorm3 /= std(unorm3)
    x3 = [unorm3 p];
    C3 = cov(x3);
    cn3 = C3[1,2:end]/norm(C3[1,2:end]);
    c3_labels = [cn3 args]

    # From the svd1/3 scatter, we see the vector that most separates clusters is svd1-svd3
    bifn_vec = cn-cn3;
    bifn_vec = bifn_vec/norm(bifn_vec);
    bifn_labels = [bifn_vec args]
 
     
#    figure()
#    p = copy(results["params"]);
#    p = p[!vec(nanrows),:];
#    paramC = cov(p);
#    vals, vecs = eig(paramC);
#    paramx = vecs[:,end-1:end]'*p';
#    scatter(paramx[1,:], paramx[2,:]);
#
#    hessians = load(farmdir*farm_id*"_hessians.jld")
#    hessians = hessians["hessians"];
#    hessians = hessians[!vec(nanrows),:,:];    
#    scalevecs = vecs[:,end-1:end]';#.*stdp;
#    scaleHess = Array(Float64, size(hessians,1),2,2);
#    for i=1:size(hessians,1)
#        scaleHess[i,:,:] = scalevecs*hessians[1,:,:]*scalevecs';
#        scaleC = inv(scaleHess[i,:,:]);
#        
#    end


 
    return F, nanrows, r_all
end

"""
    SVD_interactive(;threshold =-0.0002, plot_option=1, plot_bad_farms=true, compute_good_only=false)

Puts up an interactive plot of runs plotted in parameter SVD space. (The SVD space is defined
based on voltage traces versus time, averaged over trials.) Clicking on a dot brings up, in a
different figure, example trials from the corresponding run.

# PARAMETERS:
    
- farm_id which farm to analyze

# OPTIONAL PARAMETERS:

- threshold  training costs below this value are considered "successful" (red dots), 
             above it are "unsuccessful" (blue dots)

- plot_option  whether to use `plot_farm()` or `plot_farm2()`

- plot_bad_farms  If true, dots for the unsuccessful farms are shown, otherwise not

- compute_good_only  If true, only the successful runs are used to compute the SVD space

- opto_conditions Number of opto conditions (control only = 1)

# RETURNS

None

    
"""
function SVD_interactive(farm_id;farmdir="MiniFarms", threshold =-0.0002, plot_option=1, plot_bad_farms=true, compute_good_only=false, opto_conditions = 1)
    # get response matrix
    if opto_conditions > 1
    response, results = load(farmdir*farm_id*"_SVD_response_matrix"*string(opto_conditions)*".jld", "response","results")
    else
    response, results = load(farmdir*farm_id*"_SVD_response_matrix.jld", "response","results")
    end

    # set up filter by nan
    nanrows = any(isnan(response),2);

    # filter for good farms
    tcost = results["tcost"];
    if !compute_good_only
        tcost = tcost[!vec(nanrows),:];
        disp_cost = copy(tcost);
    end
    badcost = tcost .>= threshold;
    
    # if we are computing SVD only on the good farms, update nanrows and disp_cost
    if compute_good_only
        nanrows = nanrows | badcost;
        disp_cost = tcost[!vec(nanrows),:];
    end    

    # Filter response matrix
    r_all = response[!vec(nanrows),:];
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);
    u = copy(F[:U]); 
    u1 = u[:,1];
    u2 = u[:,3];

    # Make list of just good farms
    if compute_good_only
        # nanrows filter already selects for good farms
        u1good = u1;
        u2good = u2;
    else
        # remove bad tcost farms
        u1good = u1[!vec(badcost),:];
        u2good = u2[!vec(badcost),:];
    end

    files = results["files"];
    files = files[!vec(nanrows),:];

    function mycallback(xy, r, h, ax)
        index = find(u1 .== xy[1])
        @printf("You selected farm # %d", index[1])  
        print("\n")
        filename = files[index[1]]
        print(filename)
        print("\n")
        print(disp_cost[index[1]])
        print("\n")
        print("\n")
       
        if plot_option == 1
            plot_farm(filename)
        else
            plot_farm2(filename)
        end
        print("\n")

    end   
    pygui(true)
    BP = install_nearest_point_callback(figure(100), mycallback)
    if plot_bad_farms
        plot(u[:,1],u[:,3],"bo")
    end
    plot(u1good, u2good, "ro")
    title("SVD U columns 1 and 3")
    ylabel("SVD Dim 3")
    xlabel("SVD Dim 1")   

end


"""
    SVD_interactive2(;threshold =-0.0002, plot_option=1, plot_bad_farms=true, compute_good_only=false)

Just like `SVD)interactive()`, but puts up two synchronized side-by-side panels, showing first versus third
as well as second versus first from the U matrix of the SVD.

Thsi function puts up an interactive plot of runs plotted in parameter SVD space. (The SVD space 
is defined based on voltage traces versus time, averaged over trials.) Clicking on a dot brings up, 
in a different figure, example trials from the corresponding run.

PARAMETERS:

- farm_id which farm to analyze

# OPTIONAL PARAMETERS:

- threshold  training costs below this value are considered "successful" (red dots), 
             above it are "unsuccessful" (blue dots)

- plot_option  whether to use `plot_farm()` or `plot_farm2()`

- plot_bad_farms  If true, dots for the unsuccessful farms are shown, otherwise not

- compute_good_only  If true, only the successful runs are used to compute the SVD space

- opto_conditions Number of opto conditions (control only = 1)

# RETURNS

None

    
"""
function SVD_interactive2(farm_id;farmdir="MiniFarms", threshold =-0.0002, plot_option=1, plot_bad_farms=false, compute_good_only=false, opto_conditions=1)
    # get response matrix
    if opto_conditions > 1
    response, results = load(farmdir*farm_id*"_SVD_response_matrix"*string(opto_conditions)*".jld", "response","results")
    else
    response, results = load(farmdir*farm_id*"_SVD_response_matrix.jld", "response","results")
    end

    # set up filter by nan
    nanrows = any(isnan(response),2);

    # filter for good farms
    tcost = results["tcost"];
    if !compute_good_only
        tcost = tcost[!vec(nanrows),:];
        disp_cost = copy(tcost);
    end
    badcost = tcost .>= threshold;
    
    # if we are computing SVD only on the good farms, update nanrows and disp_cost
    if compute_good_only
        nanrows = nanrows | badcost;
        disp_cost = tcost[!vec(nanrows),:];
    end    

    # Filter response matrix
    r_all = response[!vec(nanrows),:];
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);
    u = copy(F[:U]); 
    u1 = u[:,1];
    u2 = u[:,2];
    u3 = u[:,3];

    # Make list of just good farms
    if compute_good_only
        # nanrows filter already selects for good farms
        u1good = u1;
        u2good = u2;
        u3good = u3;
    else
        # remove bad tcost farms
        u1good = u1[!vec(badcost),:];
        u2good = u2[!vec(badcost),:];
        u3good = u3[!vec(badcost),:];
    end
    files = results["files"];
    files = files[!vec(nanrows),:];

    function mycallback(xy, r, h, ax)
        if ax[:get_title]() == "SVD U columns 3 and 1"
            myX = 3; myY = 1
        elseif ax[:get_title]() == "SVD U columns 2 and 1"
            myX = 2; myY = 1
        else
            @printf("mycallback: Don't know this axis, returning\n"); return
        end
        if plot_option == 1
            index = find((u[:,myX] .== xy[1]) .& (u[:,myY] .== xy[2]))
        else
            # Alex: guess the .& operator isn't supported in Julia 0.5.2.
            # I know, I know, I should upgrade to 0.6....
            index = find((u[:,myX] .== xy[1]))
        end

        if length(index)==0; @printf("mycallback: Couldn't figure out the nearest point, returning\n"); return; end
        @printf("You selected farm # %d", index[1])  
        print("\n")
        filename = files[index[1]]
        print(filename)
        print("\n")
        print(disp_cost[index[1]])
        print("\n")
        print("\n")

        # Set all green dots to have the appropriate x and y data:
        objs = figure(100)[:findobj]()
        for o in objs
            if typeof(o) == PyCall.PyObject && contains(pystring(o), "Line2D") && o[:get_color]()=="g"
                if o[:axes][:get_title]() == "SVD U columns 3 and 1"
                    o[:set_xdata](u[index[1],3])
                    o[:set_ydata](u[index[1],1])
                    o[:set_visible](true)
                elseif o[:axes][:get_title]() == "SVD U columns 2 and 1"
                    o[:set_xdata](u[index[1],2])
                    o[:set_ydata](u[index[1],1])
                    o[:set_visible](true)
                end
            end                
        end
        
       
        if plot_option == 1
            plot_farm(filename)
        else
            plot_farm2(filename)
        end
        print("\n")

    end   
    pygui(true)
    figure(100); clf();
    BP = install_nearest_point_callback(figure(100), mycallback)
    subplot(1,2,1)
    if plot_bad_farms
        plot(u[:,3],u[:,1],"bo")
    end
    plot(u3good, u1good, "ro")
    title("SVD U columns 3 and 1")
    ylabel("SVD Dim 1")
    xlabel("SVD Dim 3")   
    plot(0, 0, "go")[1][:set_visible](false)

    subplot(1,2,2)
    if plot_bad_farms
        plot(u[:,2],u[:,1],"bo")
    end
    plot(u2good, u1good, "ro")
    title("SVD U columns 2 and 1")
    xlabel("SVD Dim 2")   
    plot(0, 0, "go")[1][:set_visible](false)
    remove_ytick_labels()
end


plot_farm_trials = 10

"""
    params = plot_farm(filename; testruns=nothing, fignum=3, overrideDict=Dict())

    Plots multiple trials from a single run of a farm.

# PARAMETERS

- filename    The filename of the .jld file containing the run, to be loaded

# OPTIONAL PARAMETERS

- testruns    Number of trials to run. Defaults to global plot_farm_trials.

- fignum      Figure to put the plot up in.

- overrideDict   A dictionary containing any model parameter values that will
              override any values loaded from the file.  For example
              `overrideDict = Dict(:sigma=>0.001)` will run with that value
              of sigma, no whater what the file said.

"""
function plot_farm(filename; testruns=nothing, fignum=3, overrideDict=Dict())

    if testruns == nothing; testruns = plot_farm_trials; end

    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

    pygui(true)
    figure(fignum); clf();
    
    pstrings = ["CONTROL", "DELAY OPTO", "CHOICE OPTO"]
    for period = 1:3
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        # :opto_strength=>0.3, 
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        # :opto_times=>["target_start-0.4" "target_start"],
        # :opto_times=>["target_start" "target_end"],
        # :post_target_period=>0.3,
        # :rule_and_delay_period=>1.2,
        # :dt=>0.005,
        ))

        # The plot_list should be the one we give it below, not whatever was in the stored parameters
        delete!(these_pars, :plot_list)

        pvax = subplot(4,3,period);   axisHeightChange(0.9, lock="t")
        pdax = subplot(4,3,period+3); axisHeightChange(0.9, lock="c"); 
        avax = subplot(4,3,period+6); axisHeightChange(0.9, lock="c")
        adax = subplot(4,3,period+9); axisHeightChange(0.9, lock="b")

        proVs, antiVs = run_ntrials(testruns, testruns; plot_list=[1:20;], plot_Us=false, 
            ax_set = Dict("pro_Vax"=>pvax, "pro_Dax"=>pdax, "anti_Vax"=>avax, "anti_Dax"=>adax),
        merge(make_dict(args, pars3, these_pars), overrideDict)...);

        hBP = length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
        hBA = length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)
        # @printf("period %d:  hBP=%.2f%%, hBA=%.2f%%\n\n", period, 100*hBP, 100*hBA)

        axes(pvax); title(@sprintf("%s  PRO hits = %.2f%%", pstrings[period], 100*hBP))
        axes(avax); title(@sprintf("ANTI hits = %.2f%%", 100*hBA))
        axes(pdax); remove_xtick_labels(); xlabel("")
        if period > 1
            remove_ytick_labels([pvax, pdax, avax, adax])
        end
        
        figure(fignum)[:canvas][:draw]()
        pause(0.001)
    end

    for a=1:length(args)
        myarg = args[a]; while length(myarg)<20; myarg=myarg*" "; end
        @printf("%s\t\t%g\n", myarg, pars3[a])
    end

    return pars3
end



"""
    plot_farm2(filename; testruns=nothing, fignum=3, overrideDict=Dict())

    Plots multiple trials from a single run of a farm. Stripped down version of plot_farm, is compatiable with Julia 0.5.2

# PARAMETERS

- filename    The filename of the .jld file containing the run, to be loaded

# OPTIONAL PARAMETERS

- testruns    Number of trials to run. Defaults to global plot_farm_trials.


- overrideDict   A dictionary containing any model parameter values that will
              override any values loaded from the file.  For example
              `overrideDict = Dict(:sigma=>0.001)` will run with that value
              of sigma, no whater what the file said.

# RETURNS
    
parameters from this farm run in <filename>
"""
function plot_farm2(filename; testruns=10, overrideDict=Dict())

    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

    pygui(true)
    
    pstrings = ["CONTROL", "DELAY OPTO", "CHOICE OPTO"]
    for period = 1:3
        f1 = period;
        f2 = period + length(pstrings);
 
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        # :opto_strength=>0.3, 
        :opto_times=>reshape(extra_pars[:opto_periods][period,:], 1, 2),
        # :opto_times=>["target_start-0.4" "target_start"],
        # :opto_times=>["target_start" "target_end"],
        # :post_target_period=>0.3,
        # :rule_and_delay_period=>1.2,
        # :dt=>0.005,
        ))

        # The plot_list should be the one we give it below, not whatever was in the stored parameters
        delete!(these_pars, :plot_list)
        proVs, antiVs = run_ntrials(testruns, testruns; plot_list=[1:10;], plot_Us=false,profig=f1,antifig=f2, merge(make_dict(args, pars3, these_pars), overrideDict)...);

        hBP = length(find(proVs[1,:]  .> proVs[4,:])) /size(proVs, 2)
        hBA = length(find(antiVs[4,:] .> antiVs[1,:]))/size(antiVs,2)
        # @printf("period %d:  hBP=%.2f%%, hBA=%.2f%%\n\n", period, 100*hBP, 100*hBA)
        figure(f1)
        title("pro "*pstrings[period]) 
        figure(f2)
        title("anti "*pstrings[period]) 
    end

    for a=1:length(args)
        myarg = args[a]; while length(myarg)<20; myarg=myarg*" "; end
        @printf("%s\t\t%g\n", myarg, pars3[a])
    end

    return pars3
end


"""
    plot_SVD_approx(rank, condition, F; opto_conditions)

    
# PARAMETERS

- rank  Which SVD dimension to plot (rank is a misnomer)

- condition Which condition (hit/error x pro/anti) to plot 

- F SVD object to analyze

# OPTIONAL PARAMETERS

- opto_conditions number of opto conditions in F
    
"""
function plot_SVD_approx(rank, condition, F; opto_conditions=1)
    if opto_conditions > 1
        error("not yet implemented, opto_conditions > 1")
    end
    figure()
    S = copy(F[:S]);
#    S[rank+1:end] = 0;    
#    low_rank = F[:U] * Diagonal(S)* F[:Vt];
    U = copy(F[:U])
    Vt = copy(F[:Vt])

    neural_dex = rank;
    if condition     == "hit-pro"
        cdex = 0;
    elseif condition == "miss-pro"
        cdex = 1;
    elseif condition == "hit-anti"
        cdex = 2;
    elseif condition == "miss-anti"
        cdex = 3;
    else
        @printf("Unknown Condition: %s \n", condition)
    end

    numts = Int(size(Vt,2)/4/4);
    numtsC= Int(size(Vt,2)/4);
    node1 = copy(Vt[neural_dex,cdex*numtsC+1:cdex*numtsC+numts]);
    node2 = copy(Vt[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts]);
    node3 = copy(Vt[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts]);
    node4 = copy(Vt[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts]);
    plot(-node1)
    plot(-node2)
    plot(-node3)
    plot(-node4)
    title(condition)
 #   return low_rank, U, S, Vt 
end


"""
    plot_psth(r_all, neural_dex; opto_conditions = 1)

Plots the average PSTH from farm run <neural_dex> in <r_all>   
assumes only one opto condition! 

# PARAMETERS

- r_all response matrix to plot from

- neural_dex which farm run to plot

# OPTIONAL PARAMETERS

- opto_conditions Number of opto conditions in r_all
    
"""
function plot_psth(r_all, neural_dex; opto_conditions=1)
    if opto_conditions > 1
        error("not yet implemented, opto_conditions > 1")
    end
    numts = Int(size(r_all,2)/4/4);
    numtsC= Int(size(r_all,2)/4);
    figure()
    cdex = 0;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 1;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 2;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])
    figure()
    cdex = 3;
    plot(r_all[neural_dex,cdex*numtsC+1:cdex*numtsC+numts])
    plot(r_all[neural_dex,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts])
    plot(r_all[neural_dex,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts])

end

"""
    build_hessian_dataset(farm_id)
    
For each farm run, compute the hessian of the test noise 

# OPTIONAL PARAMETERS:

- farmdir   Direction where farm runs are located
    
# RETURNS

Response matrix

"""
function build_hessian_dataset(farm_id; farmdir="MiniOptimized")

    
    results = load(farmdir*farm_id*"_results.jld");
    hessians = Array(Float64,length(results["cost"]),12,12);

    for i = 1:length(results["cost"])

        filename = results["files"][i];    
        ftraj3 = load(filename, "ftraj3")
        farm_hessian = ftraj3[2,end];
        hessians[i,:,:] = farm_hessian;

        # This recovers "cost3" - the training cost
#         func1 = (;params...)-> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false,seedrand=extra_pars[:seedrand], merge(merge(mypars, extra_pars), Dict(params))...)[1]       
#        farm_value, farm_grad, farm_hessian = keyword_vgh(func1, args, pars3);

        # This recovers "cost" - the test cost
#         func2 = (;params...)-> JJ(testruns, testruns; verbose=false,seedrand=extra_pars[:seedrand], merge(merge(mypars, extra_pars), Dict(params))...)[1]       
#        farm_value, farm_grad, farm_hessian = keyword_vgh(func1, args, pars3);

        # do a check to see if farm_value and tcost match
        @printf("%g %s %g\n",i, filename, results["tcost"][i])
    end

    myfilename = farmdir*farm_id*"_hessians.jld";
    save(myfilename, Dict("hessians"=>hessians))

    return hessians
end


