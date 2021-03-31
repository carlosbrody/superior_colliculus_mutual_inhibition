# This script makes Figure 7
# This script works on Julia 0.6.4

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)

# Panel 7A, Projecting solutions into 2d SVD space
SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = false)

# Panel 7B, example solutions
# 01 01 0115, 01 04 0031, 02 01 0194, 02 03 0144 

# Panel 7C, scatter plot of Anti to Pro on same side vs contra side
scatter_by_arg(results, args, "vW_PA", "dW_PA");
plot(vec([-3 3]), vec([-3 3]), "r--")
xlabel("vW (Anti -> Pro)",fontsize=16)
ylabel("dW (Anti -> Pro)",fontsize=16)

# Panel 7D, see ephys data repository

# Panel 7E, histogram of hw_pro
plot_labels = ["Weights between Pro units"]
plot_params = [results["params"][:,11]]
histo_params_2_just_one(plot_labels, plot_params, results["tcost"], results["cost"], results["files"]);


