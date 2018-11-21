##### just some notes on analysis of the C30 farm
# make_mini_farm("C32"; fromdirs="Farms_C32", todir="MiniC32")
# update_farm("C32", "MiniC32"; build_hessian=false, build_encoding=false)
# hessian = build_hessian_dataset("C32"; farmdir="MiniC32", update_only=false, get_full_filename=true)
# cluster_farms("C32"; farmdir="MiniC32", threshold = -0.0001)
# cp compute_clustering/MarinoCode_C32_MiniC32.jld MiniC32_C32_clusters.jld
# examples, cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001, num_steps=61)
# unilateral, test_farm_unilateral("C32","MiniC32";threshold=-0.0001)

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")
include("unilateral_analysis.jl")
include("cluster_farms.jl")
include("cluster_example.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
##results = load("MiniC32_C32_results.jld")
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)

# make a histogram of each parameter
HD = histo_params(f1["args"], results["params"], results["tcost"], results["cost"], results["files"]);

# scatter parameters
#cp compute_clustering/MarinoCode_C32_MiniC32.jld MiniC32_C32_clusters.jld
cluster_info    = load("MiniC32_C32_clusters.jld");
cluster_ids     = cluster_info["idx"];
scatter_by_arg(results, args, "hW_P", "hW_A"; cluster_ids = cluster_ids);

# SVD of dynamics
SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = true);
plot_PCA("C32"; farmdir="MiniC32", threshold=-0.0001, color_clusters=true,cluster_ids = cluster_ids)

# check psychometric plots
# #takesforever
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="standard")
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="binary")

# histogram for each cluster
histo_params_by_cluster(results, args, cluster_ids, target_cluster)
# some matlab functions exist in proanti/figure2, but you need to make the examples.jld in
#cluster_example.jl first 
include("cluster_example.jl")
cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001, num_steps=61, testruns=50);
#copy MiniC32_C32_examples to ~/Dropbox/proanti/figure2



#having computed the examples, lets look at some dimensionality
dims = plot_dimension_analysis(cluster_ids);
figure();
plot(sort(dims[:,7]))
plot(sort(90-dims[:,5]),"r")
xlabel("Sort index")
ylabel("Angle 2 between ref. delay space, sorted")

delay_response, target_response, results = build_epoch_response_matrix("C32"; farmdir="MiniC32");
delay_svd, target_svd, full_svd = epoch_SVD(delay_response, target_response, results; color_clusters=true, cluster_ids=cluster_ids);

dims_cluster_ids = ones(size(cluster_ids));
dims_cluster_ids[vec(dims[:,7] .>= 45)] = 2;
dims_cluster_ids2 = ones(size(cluster_ids));
dims_cluster_ids2[vec(dims[:,5] .<= 45)] = 2;

delay_svd, target_svd, full_svd = epoch_SVD(delay_response, target_response, results; color_clusters=true, cluster_ids=dims_cluster_ids);
dims = plot_dimension_analysis(dims_cluster_ids);

delay_svd, target_svd, full_svd = epoch_SVD(delay_response, target_response, results; color_clusters=true, cluster_ids=dims_cluster_ids2);
dims = plot_dimension_analysis(dims_cluster_ids2);


# having demonstrated we have two types of solutions, lets see how they map onto nodes
examples,results = load("MiniC32_C32_examples_50.jld","examples","results");
examples = get_synthetic_LR_trials(examples);
ref1, ref2 = get_reference(examples, results; seed=13,numruns=50)

SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = true, cluster_ids = dims_cluster_ids);

## unilateral analysis
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
# I have no idea what this structure means
ipsi_val = uni[:,2,1,4] - uni[:,2,2,4];
contra_val = uni[:,1,1,4] - uni[:,1,2,4];
figure()
plot(ipsi_val, contra_val, "bo")
plot(ipsi_val[unidex], contra_val[unidex], "go")
ylabel("Contralateral pro-anti")
xlabel("Ipsilateral pro-anti")
plot(vec([20]), vec([20]), "ro")
plot(vec([10]), vec([20]), "mo")

unilateral = load("MiniC32_C32_unilateral_by_strength.jld","uni_results")
numfarms = size(unilateral["uni"],1)
#ipsi/contra x pro/anti x control/delay/target/full
uni = unilateral["uni"].*100;
figure()
for i=1:11
subplot(4,3,i)
ipsi_val = uni[:,2,1,i] - uni[:,2,2,i];
contra_val = uni[:,1,1,i] - uni[:,1,2,i];
plot(ipsi_val, contra_val, "bo")
plot(ipsi_val[unidex], contra_val[unidex], "ro")
ylabel("Contralateral pro-anti")
xlabel("Ipsilateral pro-anti")
title(string(i))
ylim(-100,100)
xlim(-100,100)
end

# finding most discriminable dimension
uni_params = uni_results["params"]
params = results["params"]
non_params = params[.!unidex,:]

mu_uni = mean(uni_params,1);
mu_non = mean(non_params,1);
C_uni  = cov(uni_params);
C_non  = cov(non_params);
Sw = C_uni + C_non;
Sb = (mu_uni - mu_non)'*(mu_uni - mu_non);
invSw = inv(Sw);
vals, vecs = eig(invSw*Sb);
dvec = real(vecs[:,2])


uni_proj = zeros(size(uni_params,1),1);
for i=1:size(uni_params,1)
    uni_proj[i] = dot(uni_params[i,:],dvec);
end
non_proj = zeros(size(non_params,1),1);
for i=1:size(non_params,1)
    non_proj[i] = dot(non_params[i,:],dvec);
end

edges = -3:.5:3;
huni = fit(Histogram, vec(uni_proj),edges)
hnon = fit(Histogram, vec(non_proj),edges)
figure;
evec = collect(huni.edges[1]);
evec = (evec[2:end] - evec[1:end-1])*.5 + evec[1:end-1];
plot(evec, huni.weights./sum(huni.weights),"g")
evec2 = collect(hnon.edges[1]);
evec2 = (evec2[2:end] - evec2[1:end-1])*.5 + evec2[1:end-1];
plot(evec2, hnon.weights./sum(hnon.weights),"b")
ylabel("Prob(Solution | cluster)")
xlabel("Most Discriminable dimension")

sortrows([abs(dvec) args])
dprime = (mean(uni_proj) - mean(non_proj))./sqrt(.5*(var(uni_proj) + var(non_proj)));

