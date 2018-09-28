##### just some notes on analysis of the C30 farm
# make_mini_farm("C32"; fromdirs="Farms_C32", todir="MiniC32")
# update_farm("C32", "MiniC32"; build_hessian=false, build_encoding=false)
# hessian = build_hessian_dataset("C32"; farmdir="MiniC32", update_only=false, get_full_filename=true)
# cluster_farms("C32"; farmdir="MiniC32", threshold = -0.0001)
# cp compute_clustering/MarinoCode_C32_MiniC32.jld MiniC32_C32_clusters.jld
# examples, cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001, num_steps=61)

# still need to run
# unilateral, test_farm_unilateral("C32","MiniC32";threshold=-0.0001)

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")
include("unilateral_analysis.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
##results = load("MiniC32_C32_results.jld")
# or just load the good farms
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)

# make a histogram of each parameter
HD = histo_params(f1["args"], results["params"], results["tcost"], results["cost"], results["files"]);

# scatter parameters
#cp compute_clustering/MarinoCode_C32_MiniC32.jld MiniC32_C32_clusters.jld
cluster_info    = load("MiniC32_C32_clusters.jld");
cluster_ids     = cluster_info["idx"];
scatter_by_arg(results, args, "hW_P", "hW_A"; cluster_ids = cluster_ids);

# histogram for each cluster
histo_params_by_cluster(results, args, cluster_ids, target_cluster)
# some matlab functions exist in proanti/figure2, but you need to make the examples.jld in
#cluster_example.jl first 
include("cluster_example.jl")
cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001, num_steps=61);
#copy MiniC32_C32_examples to ~/Dropbox/proanti/figure2

# SVD of dynamics
SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = true);

# check psychometric plots
# #takesforever
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="standard")
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="binary")

plot_PCA("C32"; farmdir="MiniC32", threshold=-0.0001, color_clusters=true,cluster_ids = cluster_ids)

# unilateral analysis
unilateral = load("MiniC32_C32_unilateral.jld","uni_results")
numfarms = size(unilateral["uni"],1)
#ipsi/contra x pro/anti x control/delay/target/full
uni = unilateral["uni"].*100;

# make load_farm_unilateral_filter, and confirm results, then make histograms, etc
uni_results, uni_clusters, unidex, ipsi_unidex, contra_unidex, filters = load_farm_unilateral_filter("C32", "MiniC32"; filters=[90 70 90 70])

# save cluster definitions for making histogram plots in matlab
matwrite("/home/alex/Dropbox/proanti/figure2/uni_clusters.mat", Dict("uni_clusters"=>uni_clusters))
histo_params_by_cluster(results, args, uni_clusters,1)

# plot psychometric with new cluster definitions
plot_unilateral_psychometric("C32","MiniC32"; color_clusters="uni",uni_filters=[90 70 90 70])
plot_PCA("C32"; farmdir="MiniC32", threshold=-0.0001, color_clusters=true, cluster_ids=uni_clusters);

# Plot index plots
ipsi_val = uni[:,2,1,4] - uni[:,2,2,4];
contra_val = uni[:,1,1,4] - uni[:,1,2,4];
figure()
plot(ipsi_val, contra_val, "bo")
plot(ipsi_val[unidex], contra_val[unidex], "go")
ylabel("Contralateral pro-anti")
xlabel("Ipsilateral pro-anti")
plot(vec([20]), vec([20]), "ro")
plot(vec([10]), vec([20]), "mo")

