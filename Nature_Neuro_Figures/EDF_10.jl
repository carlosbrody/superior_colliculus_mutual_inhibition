# This script makes Extended Data Figure 10, a histogram for each parameter distribution
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

# make a histogram of each parameter
plot_args =[f1["args"][6];  f1["args"][9];  f1["args"][5];  f1["args"][11];  
            f1["args"][16]; f1["args"][3];  f1["args"][10]; f1["args"][13]; 
            f1["args"][4];  f1["args"][1];  f1["args"][2];  f1["args"][7]; 
            f1["args"][8];  f1["args"][12]; f1["args"][14]; f1["args"][15] ];
plot_labels = [ "sW (Pro)"; "sW (Anti)";     "hW (Anti)";    "hW (Pro)";
                "vW (Anti to Pro)"; "vW (Pro to anti)"; "dW (Anti to Pro)"; "dW (Pro to Anti)";
                "Noise";    "Anti Rule Input"; "Stimulus Input"; "Pro Bias"; 
                "Opto Strength"; "Target Period Input"; "Pro Rule Input"; "Constant Input"];
p = results["params"];
plot_params = [p[:,6] p[:,9] p[:,5] p[:,11]  p[:,16] p[:,3] p[:,10] p[:,13] abs.(p[:,4]) p[:,1] p[:,2] p[:,7] p[:,8] p[:,12] p[:,14] p[:,15] ];

# Plot
HD,fig = histo_params_2_internal(plot_labels, plot_params, results["tcost"], results["cost"], results["files"]);


