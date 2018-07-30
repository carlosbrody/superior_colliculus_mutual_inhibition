#just some notes on analysis of the C30 farm

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC30/farm_C30_Farms_C30_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
##results = load("MiniC30_C30_results.jld")
# or just load the good farms
results = load_farm_cost_filter("C30", "MiniC30"; threshold = -0.00025)

# make a histogram of each parameter
HD = histo_params(f1["args"], results["params"], results["tcost"], results["cost"], results["files"]);

# scatter parameters
cluster_info    = load("MiniC30_C30_clusters.jld");
cluster_ids     = cluster_info["idx"];
scatter_by_arg(results, args, "hW_P", "hW_A"; cluster_ids = cluster_ids);

# histogram for each cluster
histo_params_by_cluster(results, args, cluster_ids, target_cluster)

# SVD of dynamics
SVD_interactive("C30"; farmdir="MiniC30", threshold=-0.00025, disp_encoding=false, color_clusters = true);






