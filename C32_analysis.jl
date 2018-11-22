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
### just making figures

figure()
plot(dims[:,2].*100,dims[:,3].*100,"ko")
xlim(0,100); ylim(0,100);
ylabel("Target period V.E")
xlabel("Delay period V.E")


figure()
m1 = mean(dims[:,2]).*100;
m2 = mean(dims[:,3]).*100;
s1 = std(dims[:,2]).*100;
s2 = std(dims[:,3]).*100;

plot(m1,m2,"ko")
plot(vec([m1 m1]), vec([m2-s2 m2+s2]),"k")
plot(vec([m1-s1 m1+s1]), vec([m2 m2]),"k")
xlim(0,100); ylim(0,100);
ylabel("Target period V.E")
xlabel("Delay period V.E")

figure()
ma = mean(dims[:,1]).*100;
m1 = mean(dims[:,2]).*100;
m2 = mean(dims[:,3]).*100;
sa = std(dims[:,1]).*100;
s1 = std(dims[:,2]).*100;
s2 = std(dims[:,3]).*100;

plot(vec([1,2,3]),vec([ma, m1,m2]),"ko")
plot(vec([1,1]), vec([ma+sa, ma-sa]), "k");
plot(vec([2,2]), vec([m1+s1, m1-s1]), "k");
plot(vec([3,3]), vec([m2+s2, m2-s2]), "k");

xlim(0,4); 
ylim(0,100);
ylabel("Variance Explained %")
xticks([1,2,3],["Full, 3 PC","Delay, 2 PC","Target, 2 PC"])
ylabel("Target period V.E")
xlabel("Delay period V.E")


p = results["params"]
C = cov(p);
vals, vecs = eig(C);
totalvar = sum(vals);
figure()
plot(vec(collect(1:16)),cumsum(flipdim(vals,1))./totalvar.*100,"k")
xlabel("Parameter PC")
ylabel("Variance Explained")
ylim(0, 100)






examples,results = load("MiniC32_C32_examples_50.jld","examples","results");
dex = results["cost"] .<= -0.0001;
ex = examples[dex,:,:,:,:,:];
dt = 0.024;
tvec = vec(collect(0:60).*dt);

function get_index(ex,proanti,hitmiss,opto)
    choice_dex = zeros(size(ex,1)*size(ex,6), size(ex,5));
    rule_dex = zeros(size(ex,1)*size(ex,6), size(ex,5));

    count = 1;
    for i=1:size(ex,1)
    for j=1:size(ex,6)
        if (hitmiss & (ex[i,opto,proanti,1,end,j] > ex[i,opto,proanti,4,end,j])) | (!hitmiss & (ex[i,opto,proanti,1,end,j] < ex[i,opto,proanti,4,end,j]));
        choice_dex[count,:] = ex[i,opto,proanti,1,:,j] - ex[i,opto,proanti,4,:,j] ;
        rule_dex[count,:] = 0.5.*(ex[i,opto,proanti,1,:,j]+ex[i,opto,proanti,4,:,j])  - 0.5.*(ex[i,opto,proanti,3,:,j] +ex[i,opto,proanti,3,:,j]) ;
        count +=1;
        end
    end
    end
    choice_dex = choice_dex[1:count-1,:];
    rule_dex = rule_dex[1:count-1,:];

    return choice_dex, rule_dex;
end

choice_dex,     rule_dex    = get_index(ex,1,true,1);
choice_dexm,    rule_dexm   = get_index(ex,1,false,1);
choice_dexd,    rule_dexd   = get_index(ex,1,true,2);
choice_dexdm,   rule_dexdm  = get_index(ex,1,false,2);
choice_dext,    rule_dext   = get_index(ex,1,true,3);
choice_dextm,   rule_dextm  = get_index(ex,1,false,3);
achoice_dex,    arule_dex   = get_index(ex,2,true,1);
achoice_dexm,   arule_dexm  = get_index(ex,2,false,1);
achoice_dexd,   arule_dexd  = get_index(ex,2,true,2);
achoice_dexdm,  arule_dexdm = get_index(ex,2,false,2);
achoice_dext,   arule_dext  = get_index(ex,2,true,3);
achoice_dextm,  arule_dextm = get_index(ex,2,false,3);

# Choice index strength is roughly same for pro and anti
figure();
plot(tvec,vec(mean(choice_dex,1)),"k")
plot(tvec,vec(-mean(choice_dexm,1)),"k--")
plot(tvec,vec(mean(achoice_dex,1)),"r")
plot(tvec,vec(-mean(achoice_dexm,1)),"r--")
ylabel("Choice Index")
xlabel("Time")
ylim(-.7 , .7)
xlim(0, 1.5)

# Rule index is stronger for anti, and persists during choice period
figure();
plot(tvec,vec(mean(rule_dex,1)),"k")
plot(tvec,vec(mean(rule_dexm,1)),"k--")
plot(tvec,vec(-mean(arule_dex,1)),"r")
plot(tvec,vec(-mean(arule_dexm,1)),"r--")
ylabel("Rule Index")
xlabel("Time")
ylim(0 , .4)
xlim(0, 1.5)

# Delay inactivation selectively disrupts anti rule encoding 
figure();
plot(tvec,vec(mean(rule_dexd,1)),"k")
plot(tvec,vec(mean(rule_dexdm,1)),"k--")
plot(tvec,vec(-mean(arule_dexd,1)),"r")
plot(tvec,vec(-mean(arule_dexdm,1)),"r--")
ylabel("Rule Index")
xlabel("Time")
ylim(0 , .4)
xlim(0, 1.5)

# Choice inactivation does not disrupt anti rule encoding 
figure();
plot(tvec,vec(mean(rule_dext,1)),"k")
plot(tvec,vec(mean(rule_dextm,1)),"k--")
plot(tvec,vec(-mean(arule_dext,1)),"r")
plot(tvec,vec(-mean(arule_dextm,1)),"r--")
ylabel("Rule Index")
xlabel("Time")
ylim(0 , .4)
xlim(0, 1.5)






figure();
plot(tvec,vec(mean(choice_dexd,1)),"b")
plot(tvec,vec(mean(choice_dexdm,1)),"b--")
plot(tvec,vec(mean(choice_dext,1)),"r")
plot(tvec,vec(mean(choice_dextm,1)),"r--")
plot(tvec,vec(mean(choice_dex,1)),"k")
plot(tvec,vec(mean(choice_dexm,1)),"k--")
plot(vec([1 1]), vec([-.7 .7]),"k")
plot(vec([tvec[1] tvec[end]]), vec([0 0]),"k")
ylabel("Choice Index")
xlabel("Time")
ylim(-.7 , .7)

figure();
plot(tvec,vec(mean(rule_dexd,1)),"b")
plot(tvec,vec(mean(rule_dexdm,1)),"b--")
plot(tvec,vec(mean(rule_dext,1)),"r")
plot(tvec,vec(mean(rule_dextm,1)),"r--")
plot(tvec,vec(mean(rule_dex,1)),"k")
plot(tvec,vec(mean(rule_dexm,1)),"k--")
plot(vec([1 1]), vec([-.7 .7]),"k")
plot(vec([tvec[1] tvec[end]]), vec([0 0]),"k")
ylabel("Rule Index")
xlabel("Time")
ylim(-.7 , .7)

figure();
plot(tvec,vec(mean(choice_dexd,1)),"b")
plot(tvec,vec(mean(choice_dexdm,1)),"b--")
plot(tvec,vec(mean(choice_dext,1)),"r")
plot(tvec,vec(mean(choice_dextm,1)),"r--")
plot(tvec,vec(mean(choice_dex,1)),"k")
plot(tvec,vec(mean(choice_dexm,1)),"k--")
plot(vec([1 1]), vec([-.7 .7]),"k")
plot(vec([tvec[1] tvec[end]]), vec([0 0]),"k")
ylabel("Choice Index")
xlabel("Time")
ylim(-.7 , .7)

figure();
plot(tvec,vec(mean(rule_dexd,1)),"b")
plot(tvec,vec(mean(rule_dexdm,1)),"b--")
plot(tvec,vec(mean(rule_dext,1)),"r")
plot(tvec,vec(mean(rule_dextm,1)),"r--")
plot(tvec,vec(mean(rule_dex,1)),"k")
plot(tvec,vec(mean(rule_dexm,1)),"k--")
plot(vec([1 1]), vec([-.7 .7]),"k")
plot(vec([tvec[1] tvec[end]]), vec([0 0]),"k")
ylabel("Rule Index")
xlabel("Time")
ylim(-.7 , .7)





