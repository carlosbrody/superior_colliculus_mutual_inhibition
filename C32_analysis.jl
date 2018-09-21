#just some notes on analysis of the C30 farm

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
##results = load("MiniC32_C32_results.jld")
# or just load the good farms
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.00025)

# make a histogram of each parameter
HD = histo_params(f1["args"], results["params"], results["tcost"], results["cost"], results["files"]);

# scatter parameters
cluster_info    = load("MiniC32_C32_clusters.jld");
cluster_ids     = cluster_info["idx"];
scatter_by_arg(results, args, "hW_P", "hW_A"; cluster_ids = cluster_ids);

# histogram for each cluster
histo_params_by_cluster(results, args, cluster_ids, target_cluster)

# SVD of dynamics
SVD_interactive("C30"; farmdir="MiniC30", threshold=-0.00025, disp_encoding=false, color_clusters = true);

# check psychometric plots
# #takesforever
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="standard")
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="binarized")

# unilateral analysis
    unilateral = load("MiniC30_C30_unilateral.jld","uni_results")
    numfarms = size(unilateral["uni"],1)

#ipsi/contra x pro/anti x control/delay/target/full
    uni = unilateral["uni"].*100;

#
figure()
for i=1:numfarms
     plot(uni[i,1,1,1], uni[i,1,2,1], "rx")   
    plot(uni[i,2,1,4], uni[i,2,2,4], "bo")
end
xlabel("pro hit%")
ylabel("anti hit%")

# GO ipsi trials
figure()
for i=1:numfarms
     plot(uni[i,1,1,1], uni[i,1,2,1], "rx")   
    plot(uni[i,1,1,4], uni[i,2,2,4], "bo")
end
xlabel("pro hit%")
ylabel("anti hit%")
title("go ipsi")

# GO contra trials
figure()
for i=1:numfarms
     plot(uni[i,1,1,1], uni[i,1,2,1], "rx")   
    plot(uni[i,2,1,4], uni[i,1,2,4], "bo")
end
xlabel("pro hit%")
ylabel("anti hit%")
title("go contra")


