# This script makes Figure 6
# This script works in Julia 0.6.4

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")


# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)

# Panel C, Example solution
plot_farm_figure("MiniC32/farm_C32_Farms_C32_spock-brody01-02_0098.jld");

# Panel B
# Panel D

# Panel E, decoding
# If file does not exist, you need to make long trajectories:
# include("cluster_example.jl")
# cluster_example_trajectories("C32","MiniC32"; threshold = -0.0001, long=true)
include("C32_index_analysis.jl")


