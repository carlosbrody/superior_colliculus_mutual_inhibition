README="""
This file contains a set of functions for doing SVD based analysis of neural dynamics on the proanti model.

If you have a new farm to analyze, run update_farm() first!

============================

# BACKEND FUNCTIONS         (Need to run first, once for each farm, before any analysis)
update_farm()               runs all the following:
    results = load_farm_params()
    response= build_response_matrix()
    response= build_reduced_response_matrix()
    hessian = build_hessian_dataset()
    encoding, error_types = build_encoding_dataset();
    
# ANALYSIS FUNCTIONS        (Useful for thinking about how the model works)
SVD_interactive()           SVD analysis on dynamics, INTERACTIVE plots dynamics, and rule encoding
plot_PCA()                  plots farms in parameter PCA space, with error ellipses
plot_PCA_Redux()            compares two farms in parameter PCA space, with error ellipses
plot_SVD_cluster_approx()   Generates the average dynamics for a cluster of farms
SVD_psth_interactive()      Plots average PSTH, and SVD approximated PSTH
load_farm_cost_filter()     loads results matrix but filters by cost

# HELPER FUNCTIONS          (Used in some other analyses, but could be useful building blocks)
run_farm()                  Returns the average dynamics for this farmrun in each trial type
plot_farm()                 Plots dynamics of a farm (Stolen from Carlos in another file, not maintained)
plot_farm2()                A more stable, but less useful version of plot_farm()
plot_ellipse()              plots an ellipse
SVD_cluster_approx()        Computes average dynamics for a cluster of farms


# OLDER FUNCTIONS           (Not really maintained, keeping them around incase they are needed)
SVD_analysis()              Generates SVD and PCA analysis, non-interactive. 
plot_SVD_approx()           Generates a low rank approximation to all farm dynamics
plot_psth()                 plots the average neural response, use plot_farm() instead
SVD_interactive2()          Use SVD_interactive() instead
"""




#Need to include pro_anti.jl because it loads some packages. 
include("pro_anti.jl")
using HDF5


"""
    load_farm_params(; farm_id="C17", farmdir="MiniFarms", verbose=true, verbose_every=50)

Load final parameters for each farm run in <MiniFarms/farm_id>. Also loads the training and test cost for each farm. Saves a file <farmdir[1]>_<farm_id>_results.jld.

# OPTIONAL PARAMETERS:

- farm_id   Which farm to use

- farmdir   Direction where farm runs are located

- verbose   If true, displays some information while running

- verbose_every How often to display information

# RETURNS

A dict with fields "dirs", "files", "tcost", "cost", and "params"

"""
function load_farm_params(;farm_id="C17", farmdir="MiniOptimized", verbose=true, verbose_every=50);
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

    myfilename = farmdir[1]*"_"*farm_id*"_results.jld";
    save(myfilename, results)
    return results
end




"""
    run_farm(filename; testruns=200, overrideDict=Dict(),all_conditions=false)

Runs the farm at <filename>, <testruns> times. Then computes the average trajectory for each condition: hits/errors X pro/anti.
This function will only run the model farm for ONE duration of rule_and_delay_periods, target_periods, and post_target_periods. 
By default it will use the first duration in each of those model_params variables. 

# PARAMETERS:

- filename          The path to the farm to run

# OPTIONAL PARAMETERS:

- testruns          Number of runs to compute

- overrideDict      If you want to override some model parameters

- all_conditions    If TRUE, returns average trajectory for each opto condition as well 

# RETURNS

- a row vector that contains the average trajectory in each trial type and each opto condition

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
        :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
        :target_period=>these_pars[:target_periods][1], 
        :post_target_period=>these_pars[:post_target_periods][1], 
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
    return farm_response, numConditions
end

"""
    build_response_matrix(farm_id)
    
For each farm run, computes the average trajectory for hits/errors X pro/anti. Arranges these average trajectories in a #runs X length(hits/errors X pro/anti)*#timepoints matrix. Saves a file <farm_id>_SVD_response_matrix.jld 

If the farm has multiple durations of task events, it will only use the first one (i.e., rule_and_delay_periods...)

# PARAMETERS:

- farm_id           example "C17"

# OPTIONAL PARAMETERS:

- farmdir           Directory where farm runs are located
    
- all_conditions    If TRUE, compute response matrix for all opto conditions. If FALSE, just control trials

- update_only       If TRUE, then doesn't recompute farm runs that were already computed. Instead uses old response matrix. If TRUE, uses old_results for the list of already computed farm runs. 

- old_results       A dictionary from load_farm_params(), containing a list of files that were previously used to compute response_matrix. If update_only is false, this dictionary should be empty. 

# RETURNS

- Response matrix, where each row is the average dynamics of each node for each trial-type X opto_condition

"""
function build_response_matrix(farm_id; farmdir="MiniFarms", all_conditions = false,update_only=false, old_results=Dict())

    results = load(farmdir*"_"*farm_id*"_results.jld");
    response_matrix = [];
    numConditions = 1;
    if update_only
        if all_conditions
            numConditions = 3;
            myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix"*string(numConditions)*".jld";
        else
            myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix.jld";
        end
        old_response_matrix = load(myfilename, "response");
    end

    for i = 1:length(results["cost"])
        filename = results["files"][i];

        if update_only && (size(find(filename .== old_results["files"]),1) > 0)
            old_index = find(filename .== old_results["files"]);
            farm_response = old_response_matrix[old_index, :];
        else
            farm_response, numC = run_farm(filename;all_conditions=all_conditions);
            numConditions = numC;
        end

        if isempty(response_matrix)
            response_matrix = farm_response;
        else
            response_matrix = [response_matrix; farm_response];
        end
        @printf("%g %s %g\n",i, filename, results["tcost"][i])
    end

    if all_conditions & (numConditions > 1)
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix"*string(numConditions)*".jld";
    else
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix.jld";
    end
    save(myfilename, Dict("response"=>response_matrix, "results"=>results, "numConditions"=>numConditions ))

    return response_matrix
end

"""
    build_reduced_response_matrix(farm_id; farmdir="MiniFarms", opto_conditions = 3, time_to_truncate = 0.6)

Truncates the SVD response matrix. This is useful to equalize the time durations of different trial epochs. 
The original SVD response matrix must have already been computed. 

# PARAMETERS

- farm_id           Name of farm, example "C17"

- farmdir           Directory of farm, example "MiniFarms"

- opto_conditions   How many opto_conditions to use in the SVD response matrix, see build_response_matrix()

- time_to_truncate  Duration to remove from the start of each trial, in seconds. 

# RETURNS

Truncated response matrix, also saves the response matrix with to: <name of the original>*"_reduced.jld"

"""
function build_reduced_response_matrix(farm_id; farmdir="MiniFarms", opto_conditions = 3, time_to_truncate = 0.6)

    # load the already compute response matrix
    if opto_conditions > 1
    all_response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*".jld", "response","results")
    else
    all_response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix.jld", "response","results")
    end

    # splice out the rule parts of each trial
    # need truncation time points
    mypars, extra_pars, args, pars3 = load(results["files"][1], "mypars", "extra_pars", "args", "pars3")
    rule_and_delay_period   = mypars[:rule_and_delay_periods][1];
    target_period           = mypars[:target_periods][1];
    post_target_period      = mypars[:post_target_periods][1];
    @printf("Original rule_and_delay_period : %g\n",rule_and_delay_period)
    @printf("Original target_period         : %g\n",target_period)
    @printf("Original post_target_period    : %g\n",post_target_period)
    @printf("NEW rule_and_delay_period      : %g\n",rule_and_delay_period-time_to_truncate)

    ns = Int(floor(time_to_truncate/mypars[:dt]));
    nt = Int(floor(size(all_response,2)/opto_conditions/4/4));
    nn = Int(floor(size(all_response,2)/nt))
    response_matrix = copy(all_response[:,1:(nt-ns)*nn])

    # need to parse each row of all_response into each node/trial/opto/etc
    for i=1:nn
        sd = 1+(i-1)*(nt-ns);
        ed = i*(nt-ns);
        osd = 1+(i-1)*nt+ns;
        oed = i*nt;
        response_matrix[:,1+(i-1)*(nt-ns):i*(nt-ns)] = all_response[:,1+(i-1)*nt+ns:i*nt];
    end

    # save the reduced response matrix
     if opto_conditions > 1
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*"_reduced.jld";
    else
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix_reduced.jld";
    end   
    save(myfilename, Dict("response"=>response_matrix, "results"=>results, "numConditions"=>opto_conditions ))

    return response_matrix
end

"""
    SVD_analysis(farm_id; farmdir="MiniFarms", opto_conditions = 1, compute_good_only=false, threshold=-0.0002)

Plots a series of analyses based on the SVD of the average neural response

# PARAMETERS:

- farm_id               Which farm to analyze, example "C17"

 # OPTIONAL PARAMETERS:

- farmdir               Directory of farm runs to analyze
   
- opto_conditions       Number of opto_conditions in SVD_response_matrix

- compute_good_only     If TRUE, compute SVD only on good farms

- threshold             cutoff in cost for determining good farms

- use_reduced_SVD       if TRUE, uses the truncated SVD response matrix with a shortened rule period

# RETURNS

None

"""
function SVD_analysis(farm_id; farmdir="MiniFarms", opto_conditions = 3, compute_good_only=true, threshold=-0.0002,use_reduced_SVD=false)

    # Load responses from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
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
   
    return F, nanrows, r_all
end

"""
    SVD_interactive(farm_id;farmdir="MiniFarms", threshold =-0.0002, plot_option=1, plot_bad_farms=false, compute_good_only=true, opto_conditions = 3,disp_encoding = true, use_reduced_SVD=false, backend_mode=false, ndims=2)

Puts up an interactive plot of runs plotted in parameter SVD space. (The SVD space is defined
based on voltage traces versus time, averaged over trials.) Clicking on a dot brings up, in a
different figure, example trials from the corresponding run.

# PARAMETERS:
    
- farm_id               which farm to analyze, example "C17"

# OPTIONAL PARAMETERS:

- farmdir               Directory of which farms to analyze

- threshold             training costs below this value are considered "successful" (red dots), above it are "unsuccessful" (blue dots)

- plot_option           whether to use `plot_farm()` or `plot_farm2()` 

- plot_bad_farms        If TRUE, dots for the unsuccessful farms are shown, otherwise not

- compute_good_only     If TRUE, only the successful runs are used to compute the SVD space

- opto_conditions       Number of opto conditions (control only = 1)

- disp_encoding         If TRUE, prints rule encoding index in the console after clicking on a farm run

- use_reduced_SVD       if TRUE, uses truncated SVD response matrix with shortened rule period

- backend_mode          If TRUE, then does not plot anything. Instead it returns two variables, DATA, and string_IDs. In backend_mode, compute_good_only must be set to TRUE. If disp_encoding is also true, returns two encoding variables encoding and error_types as well.

- ndims                 Number of SVD dimensions to return in backend_mode

- psth_mode             If TRUE, then returns SVD matrix F, response mean m, and response r_all

# RETURNS

None, unless in backend mode. In backend mode:

- DATA, a npoints x ndims matrix with the SVD coordinates of each farm run

- string_IDS, a npoints vector of the file names for each farm run

If disp_encoding is also true:

- encoding, a npoints x opto_conditions x pro/anti x hit/miss x 3 encoding indicies matrix

- error_types,  a npoints x opto_conditions x pro/anti x hit/miss x 3 indicies. Instead of the indicies, shows percentage of trials with positive sign for each index. 


"""
function SVD_interactive(farm_id;farmdir="MiniFarms", threshold =-0.0002, plot_option=1, plot_bad_farms=false, compute_good_only=true, opto_conditions = 3,disp_encoding = true, use_reduced_SVD=false, backend_mode=false, ndims=2, psth_mode=false, color_clusters=false, cost_choice="cost", cluster_ids=[])

    if backend_mode & !compute_good_only
        error("Backend mode doesn't play nice with including the bad farms right now...")
    end
    if psth_mode & (!backend_mode | disp_encoding)
        error("PSTH mode requires BACKEND_MODE=true, and DISP_ENCODING=false")
    end

    # Load responses from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
    end

    if disp_encoding
    include("rule_encoding.jl")
    encoding, error_types = load(farmdir*"_"*farm_id*"_encoding.jld", "encoding","error_types")
    end
        
    if color_clusters
        if !compute_good_only
            error("Cluster coloring doesn't work with bad farms, set compute_good_only=true")
        end
        if isempty(cluster_ids)
        cluster_info    = load(farmdir*"_"*farm_id*"_clusters.jld");
        cluster_ids     = cluster_info["idx"];
        end
    end

    # set up filter by nan
    nanrows = any(isnan.(response),2);

    # filter for good farms
    cost = results[cost_choice];
    if !compute_good_only
        cost = cost[!vec(nanrows),:];
        disp_cost = copy(cost);
    end
    badcost = cost .>= threshold;
    
    # if we are computing SVD only on the good farms, update nanrows and disp_cost
    if compute_good_only
        nanrows = nanrows | badcost;
        disp_cost = cost[!vec(nanrows),:];
    end    

    # Filter response matrix
    r_all = response[!vec(nanrows),:];
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);
    u = copy(F[:U]); 
    u1 = u[:,1];
    u2 = u[:,2];
    
    # filter cluster id list
    if color_clusters
        # cluster list has been cost filtered, but not NAN filtered, so we have to filter the response matrix again to get the indices
        response_cost_filtered = response[!vec(badcost),:];
        cost_filtered_nans = any(isnan.(response_cost_filtered),2);

        cluster_ids_copy = copy(cluster_ids);
        cluster_ids_copy = cluster_ids_copy[!vec(cost_filtered_nans),:];

        if length(u1) != length(cluster_ids_copy)
            error("size mis-match between cluster ids and svd index")
        end 
    end

    # Make list of just good farms
    if compute_good_only
        # nanrows filter already selects for good farms
        u1good = u1;
        u2good = u2;
    else
        # remove bad cost farms
        u1good = u1[!vec(badcost),:];
        u2good = u2[!vec(badcost),:];
    end

    files = results["files"];
    files = files[!vec(nanrows),:];
    if disp_encoding
        encoding = encoding[!vec(nanrows),:,:,:,:];
        error_types = error_types[!vec(nanrows),:,:,:,:];
    end

    if backend_mode
        DATA = u[:,1:ndims]; 
        string_IDs = files;
        if disp_encoding
            return DATA, string_IDs, encoding, error_types
        elseif psth_mode
            return DATA, string_IDs, F, m, r_all
        else
            return DATA, string_IDs
        end
    else
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
           
            if disp_encoding
            display_encoding(encoding, error_types, index[1])
            end
            if plot_option == 1
            #    plot_farm(filename)
                plot_farm_figure(filename);
            else
                plot_farm2(filename)
            end
            print("\n")
    
        end   
        pygui(true)
        BP = install_nearest_point_callback(figure(100,figsize=(6,5.5)), mycallback)
#        xlim(-.1, .125)
#        ylim(-.125,.125)
        xticks([-.1, -.05, 0, 0.05, .1],("-0.10","-0.05","0","0.05","0.10"),fontsize=12)
        yticks([-.1, -.05, 0, 0.05, .1],("-0.10","-0.05","0","0.05","0.10"),fontsize=12)
        if plot_bad_farms
            plot(u[:,1],u[:,2],"bo")
        end

#        title("SVD U columns 1 and 2")
        ylabel("SVD Dimension 2",fontsize=16)
        xlabel("SVD Dimension 1",fontsize=16)      

        # plot each cluster in a different color
        if !color_clusters
            plot(u1good, u2good, "ko")
        else
            all_colors = "bgrcmyk"
            ids = sort(unique(cluster_ids_copy))
            for i=1:length(ids)
                if !isnan(ids[i])
                    dex = vec(cluster_ids_copy .== ids[i]);
                    plot(u1good[dex], u2good[dex], "o", color=string(all_colors[i]))
                end
            end
        end
    end

end


"""
    SVD_interactive2(;threshold =-0.0002, plot_option=1, plot_bad_farms=true, compute_good_only=false)

WARNING, this function hasn't been updated, you probably want to use SVD_interactive() unless you know what you are doing. 

Just like `SVD_interactive()`, but puts up two synchronized side-by-side panels, showing first versus third
as well as second versus first from the U matrix of the SVD.

This function puts up an interactive plot of runs plotted in parameter SVD space. (The SVD space 
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

- use_reduced_SVD   If TRUE, uses truncated SVD response matrix with shortened rule period

# RETURNS

None

    
"""
function SVD_interactive2(farm_id;farmdir="MiniFarms", threshold =-0.0002, plot_option=1, plot_bad_farms=false, compute_good_only=false, opto_conditions=1, use_reduced_SVD=false)

    # Load responses from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
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

WARNING, this function was stolen from Carlos in another file. This version is likely outdated...

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

        safe_axes(pvax); title(@sprintf("%s  PRO hits = %.2f%%", pstrings[period], 100*hBP))
        safe_axes(avax); title(@sprintf("ANTI hits = %.2f%%", 100*hBA))
        safe_axes(pdax); remove_xtick_labels(); xlabel("")
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
    build_hessian_dataset(farm_id; farmdir="MiniOptimized", force_compute_hessian=false,
                                   testruns=1600)
    
For each farm run, compute the hessian of the test noise. This function assumes the hessian has 
already been computed and is in each file, and also assumes that file farmdir*"_"*farm_id*"_results.jld"
exists, as produced by `load_farm_params()`. 

# OPTIONAL PARAMETERS:

- farmdir   Direction where farm runs are located

- force_compute_hessian    If true, will compute the Hessians instead of trying to load them from the files.

- testruns                 Only relevant if force_compute_hessian=true.  Number of Pro and Anti trials to use.
    
# RETURNS

Matrix of N runs X #params X #params hessians. Also saves that matrix to file farmdir*"_"*farm_id*"_hessians.jld"

"""
function build_hessian_dataset(farm_id; farmdir="MiniOptimized", force_compute_hessian=false,testruns=1600, update_only=false, old_results=Dict(), get_full_filename=false)

    if !isfile(farmdir*"_"*farm_id*"_results.jld")
        error(@sprintf("Don't have a results file %s to work with! Have you run load_farm_params() yet?", 
              farmdir*"_"*farm_id*"_results.jld"))
    end

    results = load(farmdir*"_"*farm_id*"_results.jld");
    nargs = size(results["params"],2)
    hessians = Array(Float64,length(results["cost"]),nargs,nargs);

    if update_only
        myfilename = farmdir*"_"*farm_id*"_hessians.jld";
        old_hessians = load(myfilename, "hessians")
    end

    for i = 1:length(results["cost"])
        filename = results["files"][i];    
        if get_full_filename
            minifilename = filename
#            filename="Farms_C30/farm_C30_"*minifilename[28:end];
             filename="Farms_C32/farm_C32_"*minifilename[28:end];
        end       
 
      if update_only && (size(find(filename .== old_results["files"]),1) > 0)
            old_index = find(filename .== old_results["files"]);
            farm_hessian = old_hessians[old_index, :,:];
            hessians[i,:,:] = farm_hessian;
      else
        if !force_compute_hessian
            fj = jldopen(filename)
            if exists(fj, "ftraj3")
                ftraj3 = load(filename, "ftraj3")
                farm_hessian = ftraj3[2,end];
            else
                warn(@sprintf("For file %s, ftraj3 not found-- can't load Hessian from it. Run me with force_compute_hessian=true to compute the Hessian.  Meanwhile, returning zeros\n", filename))
                farm_hessian = zeros(nargs, nargs)
            end
            close(fj)
            hessians[i,:,:] = farm_hessian;
        else
            mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")
            func1 = (;params...)-> JJ(testruns, testruns; verbose=false, seedrand=extra_pars[:seedrand], 
                                      merge(merge(mypars, extra_pars), Dict(params))...)[1]       
            farm_value, farm_grad, farm_hessian = keyword_vgh(func1, args, pars3);
            hessians[i,:,:] = farm_hessian;
        end
        end
        @printf("%g %s %g\n",i, filename, results["tcost"][i])
    end

    myfilename = farmdir*"_"*farm_id*"_hessians.jld";
    save(myfilename, Dict("hessians"=>hessians))

    return hessians

# THIS CODE COMPUTES THE HESSIAN FROM SCRATCH, JUST KEEPING IT HERE IN CASE IT IS NEEDED
# This recovers "cost3" - the training cost
# func1 = (;params...)-> JJ(mypars[:nPro], mypars[:nAnti]; verbose=false,seedrand=extra_pars[:seedrand], merge(merge(mypars, extra_pars), Dict(params))...)[1]       
# farm_value, farm_grad, farm_hessian = keyword_vgh(func1, args, pars3);
# This recovers "cost" - the test cost
# func2 = (;params...)->JJ(testruns, testruns; verbose=false,seedrand=extra_pars[:seedrand], merge(merge(mypars, extra_pars), Dict(params))...)[1]
# farm_value, farm_grad, farm_hessian = keyword_vgh(func1, args, pars3);
# do a check to see if farm_value and tcost match

end



"""
    plot_PCA(farm_id; farmdir="MiniOptimized", opto_conditions = 3, compute_good_only=true, threshold=-0.0002, deltaCost = 1e-4)
    
Computes PCA on parameters for each farm, and then estimates the hessian approximation to the likelihood surface consistent with a deltaCost increase in cost. This ellipse is based on the quadratic approximation to the likelihood surface. It ignores the gradient, and all higher order terms. 

# OPTIONAL PARAMETERS:

- farmdir               Directory where farm runs are located

- opto_conditions       Loads the results from the SVD response matrix computed on all opto conditions, or just the control. 
                        Note that this doesnt change the results of this file, it just changes where the results matrix is loaded from. 

- compute_good_only     If TRUE, computes PCA and ellipses only on farm runs with a cost lower than threshold

- threshold             The benchmark used to define good vs. bad runs

- deltaCost             The relative change in cost that the ellipses will plot. 

- use_reduced_SVD       If TRUE, uses truncated SVD with shortened rule period
    
# RETURNS

Nothing. Plots a figure. 


"""
function plot_PCA(farm_id; farmdir="MiniOptimized", opto_conditions = 3, compute_good_only=true, threshold=-0.0002, deltaCost = 0.0008333, use_reduced_SVD=false, color_clusters=false, cluster_ids=[], cost_choice="cost", uni_clusters = false)

    if color_clusters & !compute_good_only
        error("coloring clusters is only supported is compute_good_only=true")
    end

    # Load responses, results, and hessian from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
    end
    hessians = load(farmdir*"_"*farm_id*"_hessians.jld")
    hessians = hessians["hessians"];

    # need to filter out NaN rows 
    # Some farms in some conditions have no errors, so we have NaNs
    nanrows = any(isnan.(response),2);
    
    # Filter out farms with bad hessians
    badhess = Array(Bool,size(hessians,1),1);
    for i=1:size(hessians,1)
        evals, evecs = eig(hessians[i,:,:]);
        badhess[i] = any(evals .<=0);
    end 
    nanrows = nanrows | badhess;
 
    # Filter out farms with bad cost
    if compute_good_only
        cost = copy(results[cost_choice]);
        badcost = cost .>= threshold;
        nanrows = badcost | nanrows;
    end
    
    # get parameters
    p = copy(results["params"]);
    # filter parameters
    p = p[!vec(nanrows),:];

    if color_clusters
        # have to filter cost first, then get filter list on hessians and nans
        r2 = copy(response);
        r2 = r2[!vec(badcost),:];

        # now filter response by nan
        nanrows2 = any(isnan.(r2),2);
        cluster_ids_copy = copy(cluster_ids);
        cluster_ids_copy = cluster_ids_copy[!vec(nanrows2),:];

        # now filter by bad hessians
        hessians2 = copy(hessians);
        hessians2 = hessians2[!vec(badcost),:,:];
        hessians2 = hessians2[!vec(nanrows2),:,:];
        badhess2 = Array(Bool, size(hessians2,1),1);
        for i=1:size(hessians2,1)
            evals, evecs = eig(hessians2[i,:,:]);
            badhess2[i] = any(evals .<=0);
        end
        cluster_ids_copy = cluster_ids_copy[!vec(badhess2),:];
    end

    if color_clusters && length(cluster_ids_copy) != size(p,1)
        error("cluster labels and parameter matrix differ in size, check thresholding criteria are the same")
    end

    # gotta mean subtract and z-score
    p = p - repmat(mean(p,1), size(p,1),1);
    stdp = std(p,1);
    p = p./repmat(stdp,size(p,1),1);

    # compute scale factor for hessian
    sf = stdp;#1./stdp;
    sfsf = sf'*sf;
    
    # define covariance matrix
    paramC = cov(p);
    # do pca on parameter-covariance
    vals, vecs = eig(paramC);
    # find projections onto first and second dimensions
    paramx = vecs[:,end-1:end]'*p';

    # plot pca projections
    figure()
    fignum = gcf()[:number];
    if !color_clusters
        scatter(paramx[1,:], paramx[2,:]);
    else
        all_colors = "bgrcmyk";
        if uni_clusters
        all_colors = "gbrcmyk";
        end
        ids = sort(unique(cluster_ids_copy));
        for i=1:length(ids)
            if !isnan(ids[i])
                cdex = vec(cluster_ids_copy .== ids[i]);
                scatter(paramx[1,cdex], paramx[2,cdex], color=string(all_colors[i]))
            end
        end
    end
    
    # filter hessians
    hessians = hessians[!vec(nanrows),:,:];    
    scalevecs = vecs[:,end-1:end]';
    pcaHess = Array(Float64, size(hessians,1),2,2);

    for i=1:size(hessians,1)
        all_colors = "bgrcmyk";
        if uni_clusters
            all_colors = "gbrcmyk";
        end
        if color_clusters
        if !isnan(cluster_ids_copy[i])
            hessians[i,:,:,] = hessians[i,:,:].*sfsf;
            pcaHess[i,:,:] = scalevecs*hessians[i,:,:]*scalevecs';
            this_color = string(all_colors[Int64(cluster_ids_copy[i])])
            plot_ellipse(paramx[:,i],pcaHess[i,:,:],deltaCost,fignum;color = this_color)
        end
        end
    end
    xlabel("PCA Dim 1")
    ylabel("PCA Dim 2")

end

# center [x,y]
# covariance c (2x2)
# scale factor s
"""
    plot_ellipse(center, C,s,fignum)
    
Plots an ellipse with mean center, and major and minor axes given by the eigenvalues of C. s is a scale factor that termines the radius. Commented out code has an example on data 

# PARAMETERS:

- center            a vector for the center of the ellipse. eq [0; 0] is the origin

- C                 a matrix that defines the major and minor axes. 

- s                 a scale factor. If C is a inverse covariance matrix, then s=1.96 gives 95% confidence intervals

- fignum            Which figure to plot on. Wants the figure number, not the figure object. 
    
# RETURNS

Nothing. Plots a figure. 


"""
function plot_ellipse(center, C,s,fignum;color="r")
    
    # Make a mesh of points    
    p = 0:pi/100:(2*pi);

    # Find eigenvalues of C
    vals, vecs = eig(C);
    
    # find rotation angle of largest eigenvalue
    cosrotation = dot([1;0], vecs[:,end])/(norm([1;0])*norm(vecs[:,end]))
    rotation =( pi/2 - acos(cosrotation));
    
    # Make a rotation matrix based on angle of largest eigenvalue. 
    R = [sin(rotation) cos(rotation); -cos(rotation) sin(rotation)];

    # Here we define the radius of the ellipse along the major axis (x-axis), and minor axis (y-axis)
    xradius = (s*2/vals[end])^.5;
    yradius = (s*2/vals[1])^.5;

    # scale points by length and angle. cos(p), sin(p) make circle. 
    x = xradius*cos.(p);
    y = yradius*sin.(p); 

    # rotate x,y points back to original coordinate frame, and mean shift
    points = R*[x'; y'] + repmat(center,1,length(x));
    
    # Plot the ellipse
    figure(fignum)    
    plot(points[1,:],points[2,:],color)

    # # Test code on quadratic surface
    # # Make quadratic surface
    # H = [2 1; 1 1].*20;
    # func = (x) -> (0.5*x'*H*x)[1]
    #
    # # Make a mesh of x y points
    # xax = -21:0.1:21;
    # yax = -20:0.1:20;
    # X = repmat(xax', length(yax),1);
    # Y = repmat(yax,1,length(xax));
    # 
    # # Evaluate quadratic surface on mesh
    # J = zeros(length(xax), length(yax))
    # for  xi=1:length(xax)
    # for yi=1:length(yax)
    # J[xi,yi] = func([xax[xi], yax[yi]]);
    # end
    # end
    # 
    # # plot mesh surface, using the contour plot to demonstrate contour of deltaCost.
    # close("all")
    # figure(1)
    # deltaCost = 10
    # contour(X,Y,J',levels=[deltaCost])
    # 
    # # Plot ellipse with same deltaCost, and we match the contour. 
    # plot_ellipse([0;0],H,deltaCost,gcf()[:number])
    # 
    
end








"""
    plot_PCA_Redux(;farm_id="C17", farmdir="MiniOptimized", f2="MiniOptimized_Redux", opto_conditions = 3, compute_good_only=true, threshold=-0.0002, deltaCost = 0.0008333,plot_common_only=true,plot_connector=true, plot_f2_hessian=true)

plots the parameters for two farms, one of which has been reoptimized in PCA space. For each farm run, the ellipse representing the quadratic estimate of an increase of cost of <deltaCost>. 

# OPTIONAL PARAMETERS

- farm_id,              ex. "C17"

- farmdir,              ex. "MiniOptimized"

- f2,                   ex. "MiniOptimized_Redux", should be a reoptimized version of <farmdir>

- opto_conditions       how many opto conditions to use

- compute_good_only     if TRUE, threshold for bad farm runs

- threshold             the cutoff for defining good farms

- deltaCost             the increase in cost marked by the ellipses for each farm run

- plot_common_only      if TRUE only plots the farms that are included in both farm directory. 
                        Farms could be excluded for bad cost, or bad hessians

- plot_connector        if TRUE plots a dashed lines between the common farms in each directory.

- plot_f2_hessian       if TRUE plot the ellipse for the second directory

- use_reduced_SVD       If TRUE, uses truncated SVD response matrix with shortened rule period

"""
function plot_PCA_Redux(;farm_id="C17", farmdir="MiniOptimized", f2="MiniOptimized_Redux", opto_conditions = 3, compute_good_only=true, threshold=-0.0002, deltaCost = 0.0008333,plot_common_only=true,plot_connector=true, plot_f2_hessian=true, use_reduced_SVD=false)

    # Load responses, results, and hessian from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    response2, results2 = load(f2*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld","response","results");
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
    response2, results2 = load(f2*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld","response","results");
    end
    hessians = load(farmdir*"_"*farm_id*"_hessians.jld")
    hessians = hessians["hessians"];
    hessians2 = load(f2*"_"*farm_id*"_hessians.jld")
    hessians2 = hessians2["hessians"];

    # need to filter out NaN rows 
    # Some farms in some conditions have no errors, so we have NaNs
    nanrows = any(isnan(response),2);
    nanrows2 = any(isnan(response2),2);
    
    # Filter out farms with bad hessians
    badhess = Array(Bool,size(hessians,1),1);
    for i=1:size(hessians,1)
        evals, evecs = eig(hessians[i,:,:]);
        badhess[i] = any(evals .<=0);
    end 
    nanrows = nanrows | badhess;
    badhess2 = Array(Bool,size(hessians2,1),1);
    for i=1:size(hessians2,1)
        evals, evecs = eig(hessians2[i,:,:]);
        badhess2[i] = any(evals .<=0);
    end 
    nanrows2 = nanrows2 | badhess2;
 
    # Filter out farms with bad cost
    if compute_good_only
        tcost = copy(results["tcost"]);
        badcost = tcost .>= threshold;
        nanrows = badcost | nanrows;
        tcost2 = copy(results2["tcost"]);
        badcost2 = tcost2 .>= threshold;
        nanrows2 = badcost2 | nanrows2;
    end
    
    # gotta match up farms
    files1 = copy(results["files"]);
    files1 = files1[!vec(nanrows),:];
    for i=1:size(files1,1)
        files1[i] = files1[i][15:end]; 
    end
    files2 = copy(results2["files"]);
    files2 = files2[!vec(nanrows2),:];
    for i=1:size(files2,1)
        files2[i] = files2[i][21:end]; 
    end
    common_files = intersect(files1,files2)
    files1dex = Int[findfirst(files1,x) for x in common_files];
    files2dex = Int[findfirst(files2,x) for x in common_files];   

    # get parameters
    p1 = copy(results["params"]);
    p2 = copy(results2["params"]);
    # filter parameters
    p1 = p1[!vec(nanrows),:];
    p2 = p2[!vec(nanrows2),:];
    if plot_common_only
    p1 = p1[files1dex,:];
    p2 = p2[files2dex,:];
    end
    p=[p1; p2];

    # gotta mean subtract and z-score
    p = p - repmat(mean(p,1), size(p,1),1);
    stdp = std(p,1);
    p = p./repmat(stdp,size(p,1),1);
    p1 = p[1:size(p1,1),:];
    p2 = p[1+size(p1,1):end,:];

    # compute scale factor for hessian
    sf = stdp;#1./stdp;
    sfsf = sf'*sf;
    
    # define covariance matrix
    paramC = cov(p);
    # do pca on parameter-covariance
    vals, vecs = eig(paramC);
    # find projections onto first and second dimensions
    paramx1 = vecs[:,end-1:end]'*p1';
    paramx2= vecs[:,end-1:end]'*p2';
    # plot pca projections
    figure()
    fignum = gcf()[:number];
    if plot_common_only & plot_connector
        for i=1:size(paramx1,2)
        plot([paramx1[1,i];paramx2[1,i]],[paramx1[2,i];paramx2[2,i]],"k--");
        end
    end
    plot(paramx1[1,:], paramx1[2,:],"bo");
    plot(paramx2[1,:], paramx2[2,:],"ro");
    ylabel("PCA Dim 2")
    xlabel("PCA Dim 1")

    # filter hessians
    hessians = hessians[!vec(nanrows),:,:];    
    if plot_common_only
    hessians = hessians[files1dex,:,:];
    end
    scalevecs = vecs[:,end-1:end]';
    pcaHess = Array(Float64, size(hessians,1),2,2);
    hessians2 = hessians2[!vec(nanrows2),:,:];  
    if plot_common_only
    hessians2 = hessians2[files2dex,:,:];
    end
    pcaHess2 = Array(Float64, size(hessians2,1),2,2);

    for i=1:size(hessians,1)
        try
        hessians[i,:,:,] = hessians[i,:,:].*sfsf;
        pcaHess[i,:,:] = scalevecs*hessians[i,:,:]*scalevecs';
        plot_ellipse(paramx1[:,i],pcaHess[i,:,:],deltaCost,fignum;color="b")
        end
    end
    for i=1:size(hessians2,1)
        try
        hessians2[i,:,:,] = hessians2[i,:,:].*sfsf;
        pcaHess2[i,:,:] = scalevecs*hessians2[i,:,:]*scalevecs';
        if plot_f2_hessian
            plot_ellipse(paramx2[:,i],pcaHess2[i,:,:],deltaCost,fignum;color="r")
        end
        end
    end
    xlabel("PCA Dim 1")
    ylabel("PCA Dim 2")
    title("Blue Original, Red Redux")

end


# The includes like this don't really work inside a function, they get
# evaluated in Main context, not the inside of the function's context,
# so we have to put it out here:
include("rule_encoding.jl")

"""
    update_farm(farm_id, farmdir; build_hessian=true, build_encoding=true, update_only=false)

    Builds all necessary matrices and datasets for SVD and rule encoding analysis

# PARAMETERS

- farm_id        id of farm, eq "C17"

- farmdir        example "MiniOptimized"

# OPTIONAL PARAMETERS

- build_hessian    If true, calls `build_hessian_dataset()`

- build_encoding   If true, calls, `build_encoding_dataset()` 

- update_only       If true, only adds new farms. If false, recomputes for all farm runs

# EXAMPLE

```jldoctest
update_farm("C17","MiniOptimized")
```jldoctest

"""
function update_farm(farm_id, farmdir; build_hessian=true, build_encoding=true, update_only=true)

    # Get old results
    if update_only
        old_results = load(farmdir*"_"*farm_id*"_results.jld");
    else
        old_results = Dict();
    end
    @printf("Building Results matrix\n")
    results = load_farm_params(;farm_id=farm_id, farmdir=farmdir, verbose_every=1)

    @printf("Building Response matrix\n")
    response= build_response_matrix(farm_id; farmdir=farmdir, update_only=update_only, old_results=old_results)
    response= build_reduced_response_matrix(farm_id; farmdir=farmdir, opto_conditions=1)

    @printf("Building Response matrix all conditions\n")
    response= build_response_matrix(farm_id; farmdir=farmdir, all_conditions=true, update_only=update_only, old_results=old_results)
    response= build_reduced_response_matrix(farm_id; farmdir=farmdir, opto_conditions=3)

    if build_hessian
        @printf("Building Hessian matrix\n")
        hessian = build_hessian_dataset(farm_id; farmdir=farmdir, update_only=update_only, old_results=old_results)
    end

    if build_encoding
        @printf("Building Encoding matrix\n")
        encoding, error_types = build_encoding_dataset(farm_id; farmdir=farmdir, update_only=update_only, old_results=old_results);
    end
end

"""
    plot_SVD_cluster_approx(; farm_id="C17", farmdir="MiniOptimized", opto_conditions = 3, threshold=-0.0002, rank=2,opto_type=1, cluster_num = 1, avg_dynamics=true)

Plots the SVD approximation of the average cluster dynamics. Tries to load cluster ids from <farmdir><farm_id>_cluster_ids.jld, if that
fails, it will use a crude clustering based on the sign of U weights in SVD dim1 and dim2. 
Now plots all 4 trial types pro/anti X hit/miss

# OPTIONAL PARAMETERS
- farm_id           e.g. "C17"

- farmdir           eg "MiniOptimized"

- opto_conditions   either 1 for control only, or 3 for all opto conditions

- threshold         cutoff for defining good farms

- rank              rank of approximation. Can be any number up to number of good farms in dataset.
                    however rank > 2 will be basically the same because cluster weights will average for 0 in all dimensions above 1 or 2.

- opto_type         = 1, 2, or 3, meaning control, delay, or choice inactivations

- cluster_num       = 1,2,or 3, which cluster the farms in SVD space very crudely by svd-dim1 >/< 0, and svd-dim2 >/< 0.

- avg_dynamics      if TRUE, does approximation by average dynamics. If FALSE, uses SVD. Averaging slightly faster, same result. 

- use_reduced_SVD   If TRUE, uses truncated SVD response matrix with shortened rule period

# EXAMPLE

> plot_SVD_cluster_approx(; opto_type=1, trial_type="hit-pro", cluster_num=1)

"""
function plot_SVD_cluster_approx(; farm_id="C17", farmdir="MiniOptimized", opto_conditions = 3, threshold=-0.0002, rank=2,opto_type=1, cluster_num = 1, avg_dynamics=true, use_reduced_SVD=false)

    # Load responses, results, and hessian from all models
    if use_reduced_SVD
        reduced = "_reduced";
    else
        reduced = "";       
    end
    if opto_conditions > 1
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*string(opto_conditions)*reduced*".jld", "response","results")
    else
    response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix"*reduced*".jld", "response","results")
    end
   
    # NEED TO FILTER OUT BAD FARMS
    nanrows = any(isnan(response),2);
    tcost = results["tcost"];
    badcost = tcost .>= threshold;
    nanrows = nanrows | badcost;
    # Filter response matrix
    r_all = response[!vec(nanrows),:];
    m = mean(r_all,1);
    r_all = r_all - repmat(m, size(r_all,1),1);
    F = svdfact(r_all);
    
    # CRUDE CLUSTERING
    clusters = randn(size(r_all,1),3) .< .33;
    try
        clusters = load(farmdir*"_"*farm_id*"_cluster_ids.jld","clusters")
    catch
        utemp = copy(F[:U])
        clusters = randn(size(r_all,1),3) .< .33;
        clusters[:,1] = (utemp[:,2] .< 0) .& (utemp[:,1] .< 0);
        clusters[:,2] = (utemp[:,2] .< 0) .& (utemp[:,1] .> 0);
        clusters[:,3] = (utemp[:,2] .> 0) ;
    end

    figure()
    fignum = gcf()[:number]
    nodes = SVD_cluster_approx(clusters[:,cluster_num], rank, opto_type, "hit-pro",F,m,r_all; opto_conditions=opto_conditions,avg_dynamics=avg_dynamics)
    figure(fignum)
    subplot(2,2,1); title("hit-pro")
    plot(nodes[1,:],"b-");     plot(nodes[2,:],"r-")
    plot(nodes[3,:],"r--");     plot(nodes[4,:],"b--")
    ylim(0,1)

    nodes = SVD_cluster_approx(clusters[:,cluster_num], rank, opto_type, "miss-pro",F,m,r_all; opto_conditions=opto_conditions,avg_dynamics=avg_dynamics)
    figure(fignum)
    subplot(2,2,2); title("miss-pro")
    plot(nodes[1,:],"b-");     plot(nodes[2,:],"r-")
    plot(nodes[3,:],"r--");     plot(nodes[4,:],"b--")
    ylim(0,1)
    nodes = SVD_cluster_approx(clusters[:,cluster_num], rank, opto_type, "hit-anti",F,m,r_all; opto_conditions=opto_conditions,avg_dynamics=avg_dynamics)
    figure(fignum)
    subplot(2,2,3); title("hit-anti")
    plot(nodes[1,:],"b-");     plot(nodes[2,:],"r-")
    plot(nodes[3,:],"r--");     plot(nodes[4,:],"b--")
    ylim(0,1)
    nodes = SVD_cluster_approx(clusters[:,cluster_num], rank, opto_type, "miss-anti",F,m,r_all; opto_conditions=opto_conditions,avg_dynamics=avg_dynamics)
    figure(fignum)
    subplot(2,2,4); title("miss-anti")

    plot(nodes[1,:],"b-");     plot(nodes[2,:],"r-")
    plot(nodes[3,:],"r--");     plot(nodes[4,:],"b--")
    ylim(0,1)

    return nodes
end

"""
    Helper function for computing the SVD approximation for each cluster

# PARAMETERS

- cluster           a vector of TRUE or FALSE for inclusion of each farm into the cluster

- rank              : see doc for plot_SVD_cluster_approx()

- opto_type         : see doc for plot_SVD_cluster_approx()

- trial_type        : see doc for plot_SVD_cluster_approx()

- opto_conditions   : see doc for plot_SVD_cluster_approx()

- F                 SVD factorization for response of all farms

- m                 mean response vector that was subtracted off before SVD

- r_all             response matrix for all farms

- avg_dynamics      If TRUE, compute approximation by averaging dynamics in r_all. 
                    If FALSE, compute using SVD (same result)

"""
function SVD_cluster_approx(cluster, rank, opto_type, trial_type, F,m,r_all; opto_conditions=3, avg_dynamics=true)

    
    S = copy(F[:S]);
    U = copy(F[:U]);
    
    # Build the correct low rank in S
    S[rank+1:end] = 0;    

    # build the cluster U vector
    u = mean(U[cluster[:],:],1)
    #u[1] and u[2] should be the cluster center x,y location
  
    # Make the cluster approximation
    if avg_dynamics
        low_rank = mean(r_all[cluster[:],:],1)   
    else
        low_rank = u * Diagonal(S)* F[:Vt];
    end
    
    opto_dur = Int(length(low_rank)/opto_conditions)

    # Pull the right opto condition
    sdex = 1+opto_dur*(opto_type-1);
    edex = opto_dur*(opto_type);
    this_rank = low_rank[:,sdex:edex];
    this_mean = m[:,sdex:edex];

    # Pull the right trial type
    if trial_type     == "hit-pro"
        cdex = 0;
#  figure()
#    scatter(U[:,1],U[:,2])
#    plot(u[1],u[2], "ro")

    elseif trial_type == "miss-pro"
        cdex = 1;
    elseif trial_type == "hit-anti"
        cdex = 2;
    elseif trial_type == "miss-anti"
        cdex = 3;
    else
        @printf("Unknown Condition: %s \n", condition)
    end
    numts = Int(length(this_rank)/4/4);
    numtsC= Int(length(this_rank)/4);
    node_mean = this_mean[:,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts];

    node1 = this_rank[:,cdex*numtsC+1:cdex*numtsC+numts] + node_mean;
    node2 = this_rank[:,cdex*numtsC+1+1*numts:cdex*numtsC+2*numts] + node_mean;
    node3 = this_rank[:,cdex*numtsC+1+2*numts:cdex*numtsC+3*numts] + node_mean;
    node4 = this_rank[:,cdex*numtsC+1+3*numts:cdex*numtsC+4*numts] + node_mean;


   return [node1; node2; node3; node4]
end


"""
    plots farms in 2-D SVD space, when a farm is clicked on, plots average PSTH and SVD approximation to PSTH
    Plot Color scheme: Blue = pro, Red = anti. Solid = right nodes, dashed = left nodes. 
    plot(nodes[1,:],"b-");  plot(nodes[2,:],"r-");  plot(nodes[3,:],"r--"); plot(nodes[4,:],"b--");

# PARAMETERS

- farm_id               Name of farm to analyze, ie "C17"

# OPTIONAL PARAMETERS

- farmdir               Directory of farm runs to analyze

- threshold             Cutoff for good farm runs

- plot_bad_farms        If TRUE, will likely cause an error

- compute_good_only     If TRUE, only plot good farms

- opto_conditions       Number of opto conditions in SVD response matrix

- rank                  Rank of SVD approximation to plot, if rank > rank of SVD matrix, then max possible rank is used. 

- opto_type             = 1 (control), 2 (Delay), or 3 (Choice)

- use_reduced_SVD       If TRUE, uses truncated SVD matrix

"""
function SVD_PSTH_interactive(farm_id; farmdir="MiniOptimized", threshold=-0.0002, plot_bad_farms=false, compute_good_only=true, opto_conditions=3, rank=2, opto_type=1,use_reduced_SVD=false)
    close(100); close(101)
    remove_all_BPs()

    # Get SVD coordinates, and filenames, as well as full SVD factorization
    DATA, files, F, m, r_all = SVD_interactive(farm_id; farmdir = farmdir, threshold=threshold, plot_bad_farms=plot_bad_farms, compute_good_only=compute_good_only, opto_conditions=opto_conditions,use_reduced_SVD=use_reduced_SVD, backend_mode=true, ndims=2, disp_encoding=false, psth_mode=true);  
    if size(DATA,1) < rank
        rank = size(DATA,1);
    end

    # Set up callback function.should get filename, then call plot_PSTH_SVD_comparison
    function mycallback(xy, r, h, ax)
        index = find(DATA[:,1] .== xy[1])
        @printf("You selected farm # %d", index[1])  
        print("\n")
        filename = files[index[1]]
        print(filename)
        print("\n")
        print("\n")
        plot_PSTH_SVD_comparison(101, index[1],F,m,r_all ;rank=rank,opto_type=opto_type,opto_conditions=opto_conditions)
    end   

    # Plot SVD coordinates with callback function
    pygui(true)
    BP = install_nearest_point_callback(figure(100), mycallback)
    plot(DATA[:,1], DATA[:,2], "bo")
    title("SVD U columns 1 and 2")
    ylabel("SVD Dim 2")
    xlabel("SVD Dim 1")      
end


"""
    Plots a 4x2 set of dynamics for one farm.   First column is the SVD approximation of the PSTH dynamics. 
                                                Second column in the average PSTH dynamics.

# PARAMETERS

- fignum     Which figure to plot in

- fildex    Which farm in F, and r_all

- F         SVD factorization of r_all (from SVD_interactve())

- m         mean of r_all (from SVD_interactve())

- r_all     response matrix (from SVD_interactve())

# OPTIONAL PARAMETERS

- rank      Rank of SVD approximation

- opto_type = 1,2,3,== control, delay, choice inactivation

- opto_conditions   number of opto conditions in SVD response matrix

"""
function plot_PSTH_SVD_comparison(fignum, filedex,F,m,r_all ;rank=2,opto_type=1,opto_conditions=3)
    # Clear figure
    figure(fignum); clf()

    # set cluster to be one-hot coding of just the filename
    cluster             = Array{Bool}(size(r_all,1),1);
    cluster[:]          = false;
    cluster[filedex]    = true;

    # Then loop over trial_types
    types = ["hit-pro"; "miss-pro"; "hit-anti"; "miss-anti"];
    for i=1:4
        trial_type = types[i];
        #calls SVD_cluster_approx(average_dynamics=false) to get SVD approx
        nodes = SVD_cluster_approx(cluster, rank, opto_type, trial_type, F,m,r_all; opto_conditions=opto_conditions, avg_dynamics=false);
        subplot(4,2,1+(i-1)*2); title(trial_type*" SVD APPROX"); ylim(0,1)
        plot(nodes[1,:],"b-");  plot(nodes[2,:],"r-");  plot(nodes[3,:],"r--"); plot(nodes[4,:],"b--");
    
        #calls SVD_cluster_approx(average_dynamics=true) to get PSTH
        nodes = SVD_cluster_approx(cluster, rank, opto_type, trial_type, F,m,r_all; opto_conditions=3, avg_dynamics=true);
        subplot(4,2,2*i); title(trial_type*" PSTH"); ylim(0,1)
        plot(nodes[1,:],"b-");  plot(nodes[2,:],"r-");  plot(nodes[3,:],"r--"); plot(nodes[4,:],"b--");
 
    end
end

"""
    load_farm_cost_filter(farmid, farmdir; threshold, cost_choice)

    loads the results matrix for this farm, but filters based on cost

    OPTIONAL PARMETERS
    threshold       the cutoff for good farms
    cost_choice     base cutoff on tcost (training cost), or cost (test cost)
    
"""
function load_farm_cost_filter(farmid, farmdir; threshold=-0.0002, cost_choice="cost")
    # load results matrix
    results = load(farmdir*"_"*farmid*"_results.jld");

    # filter results matrix
    badcost = results[cost_choice] .>= threshold;
    r1 = Dict();
    r1["cost"]  = results["cost"][.!vec(badcost)];
    r1["tcost"] = results["tcost"][.!vec(badcost)];
    r1["files"] = results["files"][.!vec(badcost)];
    r1["dirs"]  = results["dirs"][.!vec(badcost)];
    r1["params"]= results["params"][.!vec(badcost),:]

    # return filtered matrix
    return r1
end 


# This is the plot-farm function that is used to generate all the plots used in the paper. 
function plot_farm_figure(filename; testruns=nothing, fignum=3, overrideDict=Dict())

    if testruns == nothing; testruns = plot_farm_trials; end

    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

    pygui(true)
    fig = figure(figsize=(12,8));  
    pstrings = ["No Inactivation", "Delay Inactivation", "Choice Inactivation"];
    p2strings = ["Pro \\ Target Left", "Anti \\ Target Left"];
    cProvall = [];
    cAntivall = [];
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
        proVs, antiVs,proVall, antiVall = run_ntrials(testruns, testruns; merge(make_dict(args, pars3, these_pars), overrideDict)...);
        tvec = collect(0:these_pars[:dt]:(these_pars[:rule_and_delay_period]+these_pars[:target_period]));
        tp = these_pars[:rule_and_delay_period];
        pvax = subplot(3,2,period+period-1);   
            title("Pro \\ Stimulus Left \n "*pstrings[period],fontsize=12)
            plot(tvec, proVall[2,:,:], color=(253/255,137/255, 57/255));
            plot(tvec, proVall[3,:,:], color=(253/255,200/255,100/255));
            plot(tvec, proVall[1,:,:], color=(0/255,  128/255,  0/255));
            plot(tvec, proVall[4,:,:], color=(100/255,200/255,100/255));
            plot(vec([tp tp]), vec([0,1]), color=(179/255, 179/255, 179/255));
            plot(vec([0 tvec[end]]), vec([0,0]), color=(179/255, 179/255, 179/255));
            ylabel("Activation",fontsize=12)
            xlabel("Time",fontsize=12)
            xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
            yticks([0, 0.5, 1],("0","0.5","1"),fontsize=12)
            ylim(-0.01,1)
            xlim(0,tvec[end])
        avax = subplot(3,2,period*2);          
            title("Anti \\ Stimulus Left \n "*pstrings[period],fontsize=12)       
            plot(tvec, antiVall[2,:,:],color=(253/255,136/255,57/255));
            plot(tvec, antiVall[3,:,:],color=(253/255,200/255,100/255));
            plot(tvec, antiVall[1,:,:],color=(0/255,128/255,0/255));
            plot(tvec, antiVall[4,:,:],color=(100/255,200/255,100/255));
            ylabel("Activation",fontsize=12)
            xlabel("Time",fontsize=12)
            plot(vec([tp tp]), vec([0,1]), color=(179/255, 179/255, 179/255))
            plot(vec([0 tvec[end]]), vec([0,0]), color=(179/255, 179/255, 179/255));
            xticks([0, 0.5, 1, 1.5],("0","0.5","1","1.5"),fontsize=12)
            yticks([0, 0.5, 1],("0","0.5","1"),fontsize=12)
            ylim(-.01,1)
            xlim(0,tvec[end])
        if period == 1
            cProvall = proVall;
            cAntivall = antiVall;
        end
    end
    fig[:subplots_adjust](hspace = .75,wspace=.75);


#    figure();
#    subplot(1,2,1);
#    pchoice = abs.(p[1,:,:] - p[4,:,:]);
#    prule = abs.(0.5.*(p[1,:,:] + p[4,:,:]) - 0.5.*(p[2,:,:] + p[3,:,:]));
#    achoice = abs.(a[1,:,:] - a[4,:,:]);
#    arule = abs.(0.5.*(a[1,:,:] + a[4,:,:]) - 0.5.*(a[2,:,:] + a[3,:,:]));
#    plot(pchoice,"k")
#    plot(achoice,"r")
#    subplot(1,2,2);
#    plot(prule,"k")
#    plot(arule,"r")
    return fig  # ,cProvall, cAntivall
end

# plots the SVD plot of all solutions, and plots the filename solution on it. used in figure 2 A. 
function get_svd_plot_figure(filename)
    svd_cord,names = SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = false, backend_mode=true);
    figure();
    plot(vec(svd_cord[:,1]), vec(svd_cord[:,2]),"ko")
    dex = find(names .== filename);
    plot(vec(svd_cord[dex,1]), vec(svd_cord[dex,2]),"ro")
end

# plots the scatter plot of vw_pa against dw_pa, and plots the filename solution on it. Used in figure 2 C
function get_parameter_plot_figure(filename)
    results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)
    f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
    args = f1["args"];
    scatter_by_arg(results, args, "vW_PA", "dW_PA");
    dex = find(results["files"] .== filename);
    plot(vec(results["params"][dex, 16]), vec(results["params"][dex, 10]),"ro")
end


