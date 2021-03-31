include("svd_cluster.jl")
using MAT

# simulates full trial inactivation
function full_trial_inactivation(farm_id, farmdir; testruns=10000)
farmfilejld = farmdir*"_"*farm_id*"_full_trial_inactivation_2.jld";
# load results
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)
new_opto = ["trial_start" "trial_end"; "trial_start" "trial_start+0.4"];
# Iterate over every farm,
# save accuracy pro/anti x 2 opto conditions
output = zeros(length(results["files"]), 2, 2);
for i=1:length(results["files"])
    println(string(i)*"/"*string(length(results["files"])))
    # get stuff for this farm
    filename = results["files"][i];
    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3");
    # run each condition
    for j=1:2
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(new_opto[j,:], 1, 2),
        :rule_and_delay_period=>these_pars[:rule_and_delay_periods][2], 
        :target_period=>these_pars[:target_periods][2], 
        :post_target_period=>these_pars[:post_target_periods][1]));
        proVs, antiVs, pfull, afull = run_ntrials(testruns, testruns; plot_list=[1:10;], plot_Us=false,
            merge(make_dict(args, pars3, these_pars), Dict())...); 
        hitsP  = 0.5*(1 + tanh.((proVs[1,:]-proVs[4,:,])/0.05));
        hitsA  = 0.5*(1 + tanh.((antiVs[4,:]-antiVs[1,:,])/0.05));
        # save each condition
        output[i,j,1] = mean(hitsP);
        output[i,j,2] = mean(hitsA);
    end 
end
# return everything for plotting
save(farmfilejld, Dict("output"=>output,"results_output"=>results))
end

# has nothing to do with clusters really, but makes example trajectories for each solution. Very useful!
function cluster_example_trajectories(farm_id, farmdir; threshold=-0.0001,testruns=50,num_steps=61)

farmfilemat = farmdir*"_"*farm_id*"_examples.mat";
farmfilejld = farmdir*"_"*farm_id*"_examples.jld";

# load results
#response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix3.jld", "response","results")
results = load_farm_cost_filter(farm_id, farmdir; threshold = threshold)


## load cluster labels
#cluster_info    = load(farmdir*"_"*farm_id*"_clusters.jld")
#cluster_ids     = cluster_info["idx"]
#ids             = sort(unique(cluster_ids));
#all_colors      = "bgrcmyk";

# Iterate over every farm,
# run 10 example trials pro and anti x 3 opto conditions
# save big matrix
examples = zeros(length(results["files"]), 3, 2, 4, num_steps, testruns);
for i=1:length(results["files"])
    println(string(i)*"/"*string(length(results["files"])))
    # get stuff for this farm
    filename = results["files"][i];
    mypars, extra_pars, args, pars3 = load(filename, "mypars", "extra_pars", "args", "pars3")

    # run each condition
    for j=1:3
        these_pars = merge(mypars, extra_pars);
        these_pars = merge(these_pars, Dict(
        :opto_times=>reshape(extra_pars[:opto_periods][j,:], 1, 2),
        :rule_and_delay_period=>these_pars[:rule_and_delay_periods][1], 
        :target_period=>these_pars[:target_periods][1], 
        :post_target_period=>these_pars[:post_target_periods][1]));

        proVs, antiVs, pfull, afull = run_ntrials(testruns, testruns; plot_list=[1:10;], plot_Us=false,
            merge(make_dict(args, pars3, these_pars), Dict())...);

        # save each condition
        examples[i,j,1,:,:,:] = pfull;
        examples[i,j,2,:,:,:] = afull;
    end
    
end

# return everything for plotting
save(farmfilejld, Dict("examples"=>examples,"results"=>results))
matwrite(farmfilemat, Dict("examples"=>examples,"results"=>results))
end

# This function computes the PCA dimensions during the rule and delay period, as well as the target period for example trials x. Then it computes the angle between those PCA dimensions and the reference PCAs 
function cluster_example_PCA(x; ref1=[], ref2=[],numruns=10,return_all=false,limited_delay=false)
    # do total PCA
    a   = copy(x);
    a   = reshape(a,2,4,61*numruns);
    da  = [a[1,:,:] a[2,:,:]];
    ca  = cov(da');
    valsa, vecsa = eig(ca);
    var_expla = [valsa[4]./sum(valsa) sum(valsa[3:4])./sum(valsa) sum(valsa[2:4])./sum(valsa) sum(valsa[:])./sum(valsa)];

    # For each solution, do PCA on the delay period activity and target period separately
    # 1 - 41, rule and delay period
    # 42 - 61, target perioda
    if limited_delay
     rd  = x[:,:,21:41,:];   
    else
    rd  = x[:,:,1:41,:];
    end
    t   = x[:,:,42:end,:];
    if limited_delay
    rd  = reshape(rd,2,4,21*numruns);
    else
    rd  = reshape(rd,2,4,41*numruns);   
    end
    d   = [rd[1,:,:] rd[2,:,:]]
    c   = cov(d');
    valsrd, vecsrd = eig(c);
    var_explrd = [valsrd[4]./sum(valsrd) sum(valsrd[3:4])./sum(valsrd) sum(valsrd[2:4])./sum(valsrd) sum(valsrd[:])./sum(valsrd)];

    t   = reshape(t,2,4,20*numruns);
    d2  = [t[1,:,:] t[2,:,:]]
    c2  = cov(d2');
    valst, vecst = eig(c2);
    var_explt    = [valst[4]./sum(valst) sum(valst[3:4])./sum(valst) sum(valst[2:4])./sum(valst) sum(valst[:])./sum(valst)];
    
    # Compute Angle between RD-pca and T-pca
    u,s,v = svd(vecsrd[:,3:4]'*vecst[:,3:4]);
    r1 = acos.(vec(s));

    if !isempty(ref1)
        # Compute Angle between RD-pca and ref1
        u,s,v       = svd(vecsrd[:,3:4]'*ref1);      
        s[s .>= 1]  = 1;
        s[s .<= -1] = -1;
        r2          = acos.(vec(s));
        # Compute Angle between T-pca and ref2
        u,s,v       = svd(vecst[:,3:4]'*ref2);       
        s[s .>= 1]  = 1;
        s[s .<= -1] = -1;
        r3          = acos.(vec(s));
    else
        r2  = [NaN NaN];
        r3  = [NaN NaN];
    end
    # Returns the variance explained during the whole trial in top 3 dims, rule_and_delay top 2 dims, target top 2 dims, angle between rd and target spaces, and angles between the reference spaces for rd and target
    if return_all
        return var_expla, var_explrd, var_explt, r1.*(90/(pi/2)), r2.*(90/(pi/2)),r3.*(90/(pi/2))
    else
        return var_expla[3], var_explrd[2], var_explt[2], r1.*(90/(pi/2)), r2.*(90/(pi/2)),r3.*(90/(pi/2))
    end
end


# This function iterates over all examples and computes the top PCA dimensions for the rule and delay period, as well as the target period, and then computes the angle between those dimensions and the reference dimensions. 
function check_dimensions(examples, results;threshold=-0.0001,ref1=[], ref2=[], numruns=10,return_all=false,limited_delay=false)
    
    # make index of which solutions to check
    dex = results["cost"] .<= threshold;

    # set up space
    if return_all
        dims = zeros(length(results["cost"]),18);
    else
        dims = zeros(length(results["cost"]),9);
    end

    # Iterate over solutions
    for i=1:length(results["cost"])
        if dex[i]
            va, v1, v2,r,r2,r3 = cluster_example_PCA(examples[i,1,:,:,:,:];ref1=ref1, ref2=ref2,numruns=numruns,return_all=return_all,limited_delay=limited_delay);
            dims[i,:] = [va v1 v2 r' r2' r3'];
        end
    end
    
    # return just the good solutions
    return dims[dex,:]
end

# This function uses a particular solution (seed) to compute a reference set of PCA vectors for the rule and delay period (ref1), and target periods (ref2). 
function get_reference(examples, results; seed = 1,numruns=10)
    # pick which farm to use as a reference
    x   = examples[seed,1,:,:,:,:];
    rd  = x[:,:,1:41,:];
    t   = x[:,:,42:end,:];

    # PCA on rule and delay period
    rd  = reshape(rd,2,4,41*numruns);
    d   = [rd[1,:,:] rd[2,:,:]]
    c   = cov(d');
    vals, vecs = eig(c);
    ref1 = vecs[:,3:4];

    # PCA on target period
    t   = reshape(t,2,4,20*numruns);
    d2  = [t[1,:,:] t[2,:,:]]
    c2  = cov(d2');
    vals2, vecs2 = eig(c2);
    ref2    = vecs2[:,3:4];
    
    # return top two dimensions for each period
    return ref1, ref2
end

# This function takes half the example trials and switches them to be left vs right trials. This allows symmetric computation of PCA spaces.
function get_synthetic_LR_trials(examples)
    new_examples = copy(examples);
    new_examples[:,:,:,1,:,51:end] = examples[:,:,:,4,:,51:end];
    new_examples[:,:,:,4,:,51:end] = examples[:,:,:,1,:,51:end];
    new_examples[:,:,:,2,:,51:end] = examples[:,:,:,3,:,51:end];
    new_examples[:,:,:,3,:,51:end] = examples[:,:,:,2,:,51:end];
    # only want to switch half
    return new_examples
end



## This is the top level function that plots dimension analysis
function plot_dimension_analysis(cluster_ids;threshold=-0.0001, testruns = 50, refseed=13,return_all=false,limited_delay=false)
    examples,results = load("MiniC32_C32_examples_50.jld","examples","results");
    examples = get_synthetic_LR_trials(examples);
    ref1, ref2 = get_reference(examples, results; seed=refseed,numruns=testruns)
    dims = check_dimensions(examples, results; threshold=threshold, ref1=ref1, ref2=ref2,numruns=testruns,return_all=return_all,limited_delay=limited_delay);

    if return_all
        return dims
    end
    figure();
    all_colors = "bgrcmyk";
    numclusters = sum(.!isnan.(unique(cluster_ids)));
    for i=1:numclusters
        cdex = vec(cluster_ids .== i);
        subplot(2,2,1)
        plot(dims[cdex,2],dims[cdex,3],"o",color=string(all_colors[i]))
        xlim(.5, 1);ylim(.5, 1)
        ylabel("Variance Explained during Target Period")
        xlabel("Variance Explained during Delay Period")

        subplot(2,2,2)
        plot(dims[cdex,4],dims[cdex,5],"o",color=string(all_colors[i]))
        xlim(0, 90);ylim(0, 90)
        plot(vec([0 90]), vec([0 90]), "k--")
        ylabel("Angle 1 between delay and target spaces")
        xlabel("Angle 2 between delay and target spaces")

        subplot(2,2,3)
        plot(dims[cdex,6],dims[cdex,7],"o",color=string(all_colors[i]))
        xlim(0, 90);ylim(0, 90)
        plot(vec([0 90]), vec([0 90]), "k--")
        ylabel("Angle 1 between reference delay space")
        xlabel("Angle 2 between reference delay space")

        subplot(2,2,4)
        plot(dims[cdex,8],dims[cdex,9],"o",color=string(all_colors[i]))
        xlim(0, 90); ylim(0, 90)
        plot(vec([0 90]), vec([0 90]), "k--")
        ylabel("Angle 1 between reference target space")
        xlabel("Angle 2 between reference target space")
    end 
    return dims
end


# This function loads the response matrix for the solution set, and splits it into epoch specific
#subsets for the delay and target periods
function build_epoch_response_matrix(farm_id; farmdir="MiniC32", all_conditions=true)
    # load full response matrix
    numconditions = 3;
    if all_conditions
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix3.jld";
    else
        myfilename = farmdir*"_"*farm_id*"_SVD_response_matrix.jld"; 
        numconditions=1;      
    end
    response, results = load(myfilename,"response","results")

    # get indexes for each time epoch
    # response is number solutions x steps_per_trial*opto_conditions*4_nodes*4trial_types (hit/missXpro/anti)
    mypars, extra_pars, args, pars3 = load(results["files"][1], "mypars", "extra_pars", "args", "pars3")
    rule_and_delay_period   = mypars[:rule_and_delay_periods][1];
    target_period           = mypars[:target_periods][1];
    post_target_period      = mypars[:post_target_periods][1];
    steps_per_trial         = Int(size(response,2)/(numconditions*4*4)); 
    delay_steps             = Int(floor(rule_and_delay_period/mypars[:dt]));
    target_steps            = steps_per_trial - delay_steps;
    numparses               = numconditions*4*4;
 
    # split matrix
    delay_response  = copy(response[:,1:delay_steps]);
    target_response = copy(response[:,delay_steps+1:steps_per_trial]);
    for i=1:numparses-1
        sd = steps_per_trial*i + 1;
        ed = sd + delay_steps -1;
        st = ed +1 ;
        et = ed + target_steps;
        delay_response  = [delay_response  copy(response[:,sd:ed])];
        target_response = [target_response copy(response[:,st:et])];
    end
   
    # return
    return delay_response, target_response, results;
end

function epoch_SVD(delay_response, target_response, results; opto_conditions = 3, threshold=-0.0001, color_clusters=false, cluster_ids=[])
    delay_svd  = generic_SVD(delay_response, results,threshold);
    target_svd = generic_SVD(target_response,results,threshold);
    full_svd   = SVD_interactive("C32";farmdir="MiniC32",threshold=threshold, disp_encoding=false, backend_mode=true);
    du = copy(delay_svd[:U]);
    tu = copy(target_svd[:U]);

 
    # plot
    ids = sort(unique(copy(cluster_ids)));
    if !color_clusters
        dcluster_ids = ones(size(du[:,1]));
        tcluster_ids = ones(size(tu[:,1]));
        fcluster_ids = ones(size(full_svd[1][:,1]));
    else
        badcost = results["cost"] .>= threshold;
        dresponse_cost = delay_response[!vec(badcost),:];
        tresponse_cost = target_response[!vec(badcost),:];
        response = load("MiniC32_C32_SVD_response_matrix3.jld","response");   
        fresponse_cost = response[!vec(badcost),:];
        dnan = any(isnan.(dresponse_cost),2);        
        tnan = any(isnan.(tresponse_cost),2); 
        fnan = any(isnan.(fresponse_cost),2); 
        dcluster_ids = copy(cluster_ids[!vec(dnan)]);
        tcluster_ids = copy(cluster_ids[!vec(tnan)]);
        fcluster_ids = copy(cluster_ids[!vec(fnan)]);
    end
    all_colors = "bgrcmyk"
    figure()

    for i=1:length(ids)
    if !isnan(ids[i])
        ddex = vec(dcluster_ids .== ids[i]);    
        tdex = vec(tcluster_ids .== ids[i]);    
        fdex = vec(fcluster_ids .== ids[i]);

        subplot(1,3,1);
        plot(du[ddex,1],du[ddex,2],"o", color=string(all_colors[i]));
        xlabel("Delay SVD Dim 1") 
        ylabel("Delay SVD Dim 2") 

        subplot(1,3,2);
        plot(tu[tdex,1],tu[tdex,2],"o", color=string(all_colors[i]));
        xlabel("Target SVD Dim 1") 
        ylabel("Target SVD Dim 2") 

        subplot(1,3,3);
        plot(full_svd[1][fdex,1], full_svd[1][fdex,2],"o", color=string(all_colors[i]))
        xlabel("Full SVD Dim 1") 
        ylabel("Full SVD Dim 2") 
    end
    end
    # return
    return delay_svd, target_svd, full_svd;
end

function generic_SVD(response,results, threshold)
    # filter bad solutions from both response matrices
    # bad = NaN and bad cost
    # compute SVD 
    nanrows = any(isnan.(response),2);
    cost    = results["cost"];
    badcost = cost .>= threshold;
    nanrows = nanrows .| badcost;
    disp_cost = cost[.!vec(nanrows),:];
    r_all   = response[.!vec(nanrows),:];
    m       = mean(r_all,1);
    r_all   = r_all - repmat(m, size(r_all,1),1);
    F       = svdfact(r_all);
    return F
end


# This analysis didn't work, leaving code for posterity
function SVD_examples(; threshold = -0.0001, cluster_ids=[])
    examples,results = load("MiniC32_C32_examples_50.jld","examples","results");
#    examples = get_synthetic_LR_trials(examples);

    # Compute full trial
    response = reshape(examples, size(examples,1),prod(size(examples)[2:end])); ;
    F = generic_SVD(response, results, threshold);
    fu = F[:U];

    dexamples = examples[:,:,:,:,1:41,:];
    texamples = examples[:,:,:,:,42:end,:];
    dresponse = reshape(dexamples, size(dexamples,1),prod(size(dexamples)[2:end])); ;
    D = generic_SVD(dresponse, results, threshold);
    du = D[:U];
    tresponse = reshape(texamples, size(texamples,1),prod(size(texamples)[2:end])); ;
    T = generic_SVD(tresponse, results, threshold);
    tu = T[:U];

    all_colors = "bgrcmyk"
    figure()
    ids = sort(unique(dims_cluster_ids));
    for i=1:length(ids)
        if !isnan(ids[i])
        dex = vec(dims_cluster_ids .== ids[i]);    

        subplot(1,3,1);
        plot(du[dex,1],du[dex,2],"o", color=string(all_colors[i]));
        xlabel("Delay Example SVD Dim 1") 
        ylabel("Delay Example SVD Dim 2") 

        subplot(1,3,2);
        plot(tu[dex,1],tu[dex,2],"o", color=string(all_colors[i]));
        xlabel("Target Example SVD Dim 1") 
        ylabel("Target Example SVD Dim 2") 

        subplot(1,3,3);
        plot(fu[dex,1], fu[dex,2],"o", color=string(all_colors[i]))
        xlabel("Full Example SVD Dim 1") 
        ylabel("Full Example SVD Dim 2") 

        end
    end

    ids = sort(unique(dims_cluster_ids));
    for i=1:length(ids)
        if !isnan(ids[i])
        dex = vec(dims_cluster_ids .== ids[i]);    

        figure(1);
        plot3D(du[dex,1],du[dex,2],du[dex,3],"o", color=string(all_colors[i]));
        figure(2);
        plot3D(tu[dex,1],tu[dex,2],tu[dex,3],"o", color=string(all_colors[i]));
        figure(3);
        plot3D(fu[dex,1],fu[dex,2],fu[dex,3],"o", color=string(all_colors[i]));
        end
    end
end



