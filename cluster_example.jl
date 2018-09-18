include("svd_cluster.jl")
using MAT

function cluster_example_trajectories(farm_id, farmdir; threshold=-0.00025,testruns=10,num_steps=76)

farmfilemat = farmdir*"_"*farm_id*"_examples.mat";
farmfilejld = farmdir*"_"*farm_id*"_examples.jld";

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



