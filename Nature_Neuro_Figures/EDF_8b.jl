# This script makes the EDF 8b
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

# variance explained across parameter dimensions and dynamics
p = results["params"]
C = cov(p);
vals, vecs = eig(C);
totalvar = sum(vals);
output = SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = false, backend_mode=true, psth_mode=true);
F = output[3];
s = F[:S];
ve_SVD = cumsum(s.^2./sum(s.^2)).*100;

figure(figsize=(3.5,3))
plot(vec(collect(1:16)),cumsum(flipdim(vals,1))./totalvar.*100,"k")
plot(vec(collect(1:16)),vec(ve_SVD[1:16]),"r")
xlabel("# Dimensions",fontsize=16)
ylabel("Variance Explained\nAcross Solutions",fontsize=16)
yticks(fontsize=12)
xticks([5,10,15],("5","10","15"),fontsize=12)
xlim(1,16)
ylim(0, 100)
plt[:tight_layout]()

