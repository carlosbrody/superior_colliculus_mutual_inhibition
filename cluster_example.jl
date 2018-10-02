include("svd_cluster.jl")
using MAT

function cluster_example_trajectories(farm_id, farmdir; threshold=-0.00025,testruns=10,num_steps=76)

farmfilemat = farmdir*"_"*farm_id*"_examples_50.mat";
farmfilejld = farmdir*"_"*farm_id*"_examples_50.jld";

# load results
response, results = load(farmdir*"_"*farm_id*"_SVD_response_matrix3.jld", "response","results")

# load cluster labels
cluster_info    = load(farmdir*"_"*farm_id*"_clusters.jld")
cluster_ids     = cluster_info["idx"]
ids             = sort(unique(cluster_ids));
all_colors      = "bgrcmyk";

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
save(farmfilejld, Dict("examples"=>examples,"results"=>results, "cluster_ids"=>cluster_ids))
matwrite(farmfilemat, Dict("examples"=>examples,"results"=>results, "cluster_ids"=>cluster_ids))
end


function cluster_example_PCA(x; ref1=[], ref2=[],numruns=10)
    # do total PCA
    a = copy(x);
    a = reshape(a,2,4,61*numruns);
    da = [a[1,:,:] a[2,:,:]];
    ca = cov(da');
    valsa, vecsa = eig(ca);
    var_expla = [valsa[4]./sum(valsa) sum(valsa[3:4])./sum(valsa) sum(valsa[2:4])./sum(valsa) sum(valsa[:])./sum(valsa)];
    # For each solution, do PCA on the delay period activity and target period separately
    # 1 - 41, rule and delay period
    # 42 - 61, target perioda
    rd = x[:,:,1:41,:];
    t  = x[:,:,42:end,:];
    rd= reshape(rd,2,4,41*numruns);
    d = [rd[1,:,:] rd[2,:,:]]
    c = cov(d');
    vals, vecs = eig(c);
    var_expl = [vals[4]./sum(vals) sum(vals[3:4])./sum(vals) sum(vals[2:4])./sum(vals) sum(vals[:])./sum(vals)];

    t= reshape(t,2,4,20*numruns);
    d2 = [t[1,:,:] t[2,:,:]]
    c2 = cov(d2');
    vals2, vecs2 = eig(c2);
    var_expl2 = [vals2[4]./sum(vals2) sum(vals2[3:4])./sum(vals2) sum(vals2[2:4])./sum(vals2) sum(vals2[:])./sum(vals2)];
    
    u,s,v = svd(vecs[:,3:4]'*vecs2[:,3:4]);

    if !isempty(ref1)
        u,s,v = svd(vecs[:,3:4]'*ref1);      
        s[s .>= 1] = 1;
        s[s .<= -1] = -1;
        r2 = acos.(vec(s));
        u,s,v = svd(vecs2[:,3:4]'*ref2);       
        s[s .>= 1] = 1;
        s[s .<= -1] = -1;
        r3 = acos.(vec(s));
    else
    r2 = [];
    r3=[];
    end
    # Returns the variance explained during the whole trial, rule_and_delay, target, angle between rd and target spaces, and angles between the reference spaces for rd and target
    return var_expla[3], var_expl[2], var_expl2[2], acos.(vec(s)).*(90/(pi/2)), r2.*(90/(pi/2)),r3.*(90/(pi/2))
end

function check_dimensions(examples, results;threshold=-0.0001,ref1=[], ref2=[], numruns=10)
    
    # make index of which solutions to check
    dex = results["cost"] .<= threshold;

    # set up space
    dims = zeros(length(results["cost"]),9);

    # Iterate over solutions
    for i=1:length(results["cost"])
        if dex[i]
            va, v1, v2,r,r2,r3 = cluster_example_PCA(examples[i,1,:,:,:,:];ref1=ref1, ref2=ref2,numruns=numruns);
            dims[i,:] = [v1 v2 r' r2' r3' va];
        end
    end
    
    # return just the good solutions
    return dims[dex,:]
end

function get_reference(examples, results; seed = 1,numruns=10)
    # pick which farm to use as a reference
    x = examples[seed,1,:,:,:,:];

    # PCA on rule and delay period
    rd = x[:,:,1:41,:];
    t  = x[:,:,42:end,:];
    rd= reshape(rd,2,4,41*numruns);
    d = [rd[1,:,:] rd[2,:,:]]
    c = cov(d');
    vals, vecs = eig(c);
    ref1 = vecs[:,3:4];

    # PCA on target period
    t= reshape(t,2,4,20*numruns);
    d2 = [t[1,:,:] t[2,:,:]]
    c2 = cov(d2');
    vals2, vecs2 = eig(c2);
    ref2 = vecs2[:,3:4];
    
    # return top two dimensions for each period
    return ref1, ref2
end

function get_synthetic_LR_trials(examples)
    new_examples = copy(examples);
    new_examples[:,:,:,1,:,26:end] = examples[:,:,:,4,:,26:end];
    new_examples[:,:,:,4,:,26:end] = examples[:,:,:,1,:,26:end];
    new_examples[:,:,:,2,:,26:end] = examples[:,:,:,3,:,26:end];
    new_examples[:,:,:,3,:,26:end] = examples[:,:,:,2,:,26:end];

    # only want to switch half
    return new_examples
end









