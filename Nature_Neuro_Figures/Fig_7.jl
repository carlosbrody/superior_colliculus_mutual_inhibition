# keep notes on which examples are used in the paper figures
#    Fig 1   01-02_0098
#    Fig 2   01 01 0115, 01 04 0031, 02 01 0194, 02 03 0144 

# load relevant files
include("results_analysis.jl")
include("svd_cluster.jl")
include("parameter_analysis.jl")

# load (arbitrary first farm for list of parameter names "args")
f1 = load("MiniC32/farm_C32_Farms_C32_spock-brody01-01_0001.jld")
args = f1["args"];

# load results of entire farm run
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)

# Panel 7C, scatter plot of Anti to Pro on same side vs contra side
scatter_by_arg(results, args, "vW_PA", "dW_PA");
plot(vec([-3 3]), vec([-3 3]), "r--")
xlabel("vW (Anti -> Pro)",fontsize=16)
ylabel("dW (Anti -> Pro)",fontsize=16)




# SVD of dynamics. This makes the figure used in Model figure 2 panel A.
SVD_interactive("C32"; farmdir="MiniC32", threshold=-0.0001, disp_encoding=false, color_clusters = false);
plot_PCA("C32"; farmdir="MiniC32", threshold=-0.0001, color_clusters=true,cluster_ids = cluster_ids)

# check psychometric plots. This plots the performance on a test set of noise
plot_psychometric(results, color_clusters=true, cluster_ids=cluster_ids, plot_only=Inf, hit_type="standard")


# histogram for each cluster
histo_params_by_cluster(results, args, cluster_ids, target_cluster)
# two better matlab functions exist for doing this, but requires two steps
#1). cluster_example_trajectories("C32", "MiniC32"; threshold=-0.0001, num_steps=61, testruns=50);
#2) copy MiniC32_C32_examples to ~/Dropbox/proanti/figure2
#3) run the matlab files


#having computed the examples, lets look at some dimensionality
# This code looks at the dimensionality of each solutions dynamics, and looks at angles between relevant subspaces. 
# interesting idea, we never used it, probably needs an updating. 
dims = plot_dimension_analysis(cluster_ids);
figure();
plot(sort(dims[:,7]))
plot(sort(90-dims[:,5]),"r")
xlabel("Sort index")
ylabel("Angle 2 between ref. delay space, sorted")
delay_response, target_response, results = build_epoch_response_matrix("C32"; farmdir="MiniC32");
delay_svd, target_svd, full_svd = epoch_SVD(delay_response, target_response, results; color_clusters=true, cluster_ids=cluster_ids);

##### All the unilateral analysis should be treated with extreme caution. Alex believes there was a bug somewhere. 
## unilateral analysis
# unilateral = load("MiniC32_C32_unilateral.jld","uni_results")
# numfarms = size(unilateral["uni"],1)
# #ipsi/contra x pro/anti x control/delay/target/full
# uni = unilateral["uni"].*100;
# 
# # make load_farm_unilateral_filter, and confirm results, then make histograms, etc
# uni_results, uni_clusters, unidex, ipsi_unidex, contra_unidex, filters = load_farm_unilateral_filter("C32", "MiniC32"; filters=[90 70 90 70])
# 
# # save cluster definitions for making histogram plots in matlab
# matwrite("/home/alex/Dropbox/proanti/figure2/uni_clusters.mat", Dict("uni_clusters"=>uni_clusters))
# histo_params_by_cluster(results, args, uni_clusters,1)
# 
# # plot psychometric with new cluster definitions
# plot_unilateral_psychometric("C32","MiniC32"; color_clusters="uni",uni_filters=[90 70 90 70])
# plot_PCA("C32"; farmdir="MiniC32", threshold=-0.0001, color_clusters=true, cluster_ids=uni_clusters);
# Plot index plots
# I have no idea what this structure means
# ipsi_val = uni[:,2,1,4] - uni[:,2,2,4];
# contra_val = uni[:,1,1,4] - uni[:,1,2,4];
# figure()
# plot(ipsi_val, contra_val, "bo")
# plot(ipsi_val[unidex], contra_val[unidex], "go")
# ylabel("Contralateral pro-anti")
# xlabel("Ipsilateral pro-anti")
# plot(vec([20]), vec([20]), "ro")
# plot(vec([10]), vec([20]), "mo")
# 
# unilateral = load("MiniC32_C32_unilateral_by_strength.jld","uni_results")
# numfarms = size(unilateral["uni"],1)
# #ipsi/contra x pro/anti x control/delay/target/full
# uni = unilateral["uni"].*100;
# figure()
# for i=1:11
# subplot(4,3,i)
# ipsi_val = uni[:,2,1,i] - uni[:,2,2,i];
# contra_val = uni[:,1,1,i] - uni[:,1,2,i];
# plot(ipsi_val, contra_val, "bo")
# plot(ipsi_val[unidex], contra_val[unidex], "ro")
# ylabel("Contralateral pro-anti")
# xlabel("Ipsilateral pro-anti")
# title(string(i))
# ylim(-100,100)
# xlim(-100,100)
# end
# # finding most discriminable dimension
# uni_params = uni_results["params"]
# params = results["params"]
# non_params = params[.!unidex,:]
# mu_uni = mean(uni_params,1);
# mu_non = mean(non_params,1);
# C_uni  = cov(uni_params);
# C_non  = cov(non_params);
# Sw = C_uni + C_non;
# Sb = (mu_uni - mu_non)'*(mu_uni - mu_non);
# invSw = inv(Sw);
# vals, vecs = eig(invSw*Sb);
# dvec = real(vecs[:,2])
# uni_proj = zeros(size(uni_params,1),1);
# for i=1:size(uni_params,1)
#     uni_proj[i] = dot(uni_params[i,:],dvec);
# end
# non_proj = zeros(size(non_params,1),1);
# for i=1:size(non_params,1)
#     non_proj[i] = dot(non_params[i,:],dvec);
# end
# 
# edges = -3:.5:3;
# huni = fit(Histogram, vec(uni_proj),edges)
# hnon = fit(Histogram, vec(non_proj),edges)
# figure;
# evec = collect(huni.edges[1]);
# evec = (evec[2:end] - evec[1:end-1])*.5 + evec[1:end-1];
# plot(evec, huni.weights./sum(huni.weights),"g")
# evec2 = collect(hnon.edges[1]);
# evec2 = (evec2[2:end] - evec2[1:end-1])*.5 + evec2[1:end-1];
# plot(evec2, hnon.weights./sum(hnon.weights),"b")
# ylabel("Prob(Solution | cluster)")
# xlabel("Most Discriminable dimension")
# 
# sortrows([abs(dvec) args])
# dprime = (mean(uni_proj) - mean(non_proj))./sqrt(.5*(var(uni_proj) + var(non_proj)));
### just making figures

#####
#####
# making psychometrics figures for model figure 1, panels B and D
#####
#####
results = load_farm_cost_filter("C32", "MiniC32"; threshold = -0.0001)
output = plot_psychometric(results, plot_only=Inf, hit_type="standard",backend_mode=true)
full_output = load("MiniC32_C32_full_trial_inactivation.jld","output");
full_output = load("MiniC32_C32_full_trial_inactivation_2.jld","output");
# Making psychometric with relative error increase
## Delay period
targets = [90 70; 85 50; 90 70];
diff_targets = [5 20; 0 0];
mean_diff_pro_delay = mean(output[:,2,1] - output[:,1,1]).*-100;
mean_diff_anti_delay = mean(output[:,2,2] - output[:,1,2]).*-100;
std_diff_pro_delay = std(output[:,2,1] - output[:,1,1]).*-100;
std_diff_anti_delay = std(output[:,2,2] - output[:,1,2]).*-100;
figure(figsize=(1.5,2))
plot(vec([1 1.5]), vec([diff_targets[1,1] diff_targets[1,1]]),"k")
plot(vec([2 2.5]), vec([diff_targets[1,2] diff_targets[1,2]]),"k")
plot(vec([1.25]), mean_diff_pro_delay, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_delay, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_delay + 1.96*[-std_diff_pro_delay, std_diff_pro_delay]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_delay + 1.96*[-std_diff_anti_delay, std_diff_anti_delay]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-5, 25)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Delay")
xticks([1.25, 2.25], ["Pro", "Anti"])
## Choice period
mean_diff_pro_choice = mean(output[:,3,1] - output[:,1,1]).*-100;
mean_diff_anti_choice = mean(output[:,3,2] - output[:,1,2]).*-100;
std_diff_pro_choice = std(output[:,3,1] - output[:,1,1]).*-100;
std_diff_anti_choice = std(output[:,3,2] - output[:,1,2]).*-100;
figure(figsize=(1.5,2))
plot(vec([1 1.5]), vec([diff_targets[2,1] diff_targets[2,1]]),"k")
plot(vec([2 2.5]), vec([diff_targets[2,2] diff_targets[2,2]]),"k")
plot(vec([1.25]), mean_diff_pro_choice, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_choice, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_choice + 1.96*[-std_diff_pro_choice, std_diff_pro_choice]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_choice + 1.96*[-std_diff_anti_choice, std_diff_anti_choice]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-5, 25)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Choice")
xticks([1.25, 2.25], ["Pro", "Anti"])
## Full Trial
if true
diff_pro_full  = (full_output[:,1,1]    - .9).*-100;
diff_anti_full = (full_output[:,1,2]   - .7).*-100;
plb = mean(diff_pro_full) - 3*std(diff_pro_full);
pub = mean(diff_pro_full) + 3*std(diff_pro_full);
pdex = (diff_pro_full .> plb) .& (diff_pro_full .< pub);
alb = mean(diff_anti_full) - 3*std(diff_anti_full);
aub = mean(diff_anti_full) + 3*std(diff_anti_full);
adex = (diff_anti_full .> alb) .& (diff_anti_full .< aub);
dex = adex .* pdex;
mean_diff_pro_full = mean(full_output[dex,1,1]    - .9).*-100;
mean_diff_anti_full = mean(full_output[dex,1,2]   - .7).*-100;
std_diff_pro_full = std(full_output[dex,1,1]      - .9).*-100;
std_diff_anti_full = std(full_output[dex,1,2]     - .7).*-100;
else
mean_diff_pro_full = mean(full_output[:,1,1]    - .9).*-100;
mean_diff_anti_full = mean(full_output[:,1,2]   - .7).*-100;
std_diff_pro_full = std(full_output[:,1,1]      - .9).*-100;
std_diff_anti_full = std(full_output[:,1,2]     - .7).*-100;
end
figure(figsize=(1.5,2))
plot(vec([1.25]), mean_diff_pro_full, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_full, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_full  + 1.96*[-std_diff_pro_full,  std_diff_pro_full]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_full + 1.96*[-std_diff_anti_full, std_diff_anti_full]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-15, 35)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Full trial")
xticks([1.25, 2.25], ["Pro", "Anti"])
## Rule period
if true
diff_pro_rule  = (full_output[:,2,1]    - .9).*-100;
diff_anti_rule = (full_output[:,2,2]   - .7).*-100;
plb = mean(diff_pro_rule) - 3*std(diff_pro_rule);
pub = mean(diff_pro_rule) + 3*std(diff_pro_rule);
pdex = (diff_pro_rule .> plb) .& (diff_pro_rule .< pub);
alb = mean(diff_anti_rule) - 3*std(diff_anti_rule);
aub = mean(diff_anti_rule) + 3*std(diff_anti_rule);
adex = (diff_anti_rule .> alb) .& (diff_anti_rule .< aub);
dex = adex .* pdex;
mean_diff_pro_rule = mean(full_output[dex,2,1]    - .9).*-100;
mean_diff_anti_rule = mean(full_output[dex,2,2]   - .7).*-100;
std_diff_pro_rule = std(full_output[dex,2,1]      - .9).*-100;
std_diff_anti_rule = std(full_output[dex,2,2]     - .7).*-100;
else
mean_diff_pro_rule = mean(full_output[:,2,1]    - output[:,1,1]).*-100;
mean_diff_anti_rule = mean(full_output[:,2,2]   - output[:,1,2]).*-100;
std_diff_pro_rule = std(full_output[:,2,1]      - output[:,1,1]).*-100;
std_diff_anti_rule = std(full_output[:,2,2]     - output[:,1,2]).*-100;
end

figure(figsize=(1.5,2))
plot(vec([1.25]), mean_diff_pro_rule, "o",color=(0, 128/255, 0))
plot(vec([2.25]), mean_diff_anti_rule, "o", color=(253/255,137/255,57/255))
plot(vec([1.25 1.25]), vec([mean_diff_pro_rule  + 1.96*[-std_diff_pro_rule,  std_diff_pro_rule]][1]),color=(0, 128/255, 0))
plot(vec([2.25 2.25]), vec([mean_diff_anti_rule + 1.96*[-std_diff_anti_rule, std_diff_anti_rule]][1]), color=(253/255,137/255,57/255))
plot(vec([0 3]), vec([0 0]), color=(179/255,179/255,179/255))
ylim(-15, 35)
xlim(0.5, 3)
ylabel("% Error Increase")
xlabel("Task cue")
xticks([1.25, 2.25], ["Pro", "Anti"])


#####
#####
# VARIANCE EXPLAINED PLOTS. These make supplementary model figure 1
#####
#####
# bar plot with full trial in 3 dims, each epoch in 2 dims
dims = plot_dimension_analysis(cluster_ids;return_all=true,limited_delay=true);
figure(figsize=(3,3))
means = mean(dims,1)[1:12].*100;
stds  = std(dims,1)[1:12].*100;
plot(vec([0.5 3.5]), vec([86 86]), color=(179/255,179/255,179/255));
plot(vec([0.5 3.5]), vec([97 97]), color=(179/255,179/255,179/255));
plot(vec([1 2 3]), vec(means[1:4:end]),"ko-")
plot(vec([1 2 3]), vec(means[2:4:end]),"ro-")
plot(vec([1 2 3]), vec(means[3:4:end]),"bo-")

plot(vec([1,1]), vec([means[1]+stds[1], means[1]-stds[1]]), "k");
plot(vec([1,1]), vec([means[2]+stds[2], means[2]-stds[2]]), "r");
plot(vec([1,1]), vec([means[3]+stds[3], means[3]-stds[3]]), "b");
plot(vec([2,2]), vec([means[5]+stds[5], means[5]-stds[5]]), "k");
plot(vec([2,2]), vec([means[6]+stds[6], means[6]-stds[6]]), "r");
plot(vec([2,2]), vec([means[7]+stds[7], means[7]-stds[7]]), "b");
plot(vec([3,3]), vec([means[9]+stds[9], means[9]-stds[9]]), "k");
plot(vec([3,3]), vec([means[10]+stds[10], means[10]-stds[10]]), "r");
plot(vec([3,3]), vec([means[11]+stds[11], means[11]-stds[11]]), "b");

xlim(0.5,3.5); 
ylim(0,100);
ylabel("Variance Explained %",fontsize=16)
yticks(fontsize=12)
xticks([1,2,3],["Full","Delay","Target"],fontsize=12)


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

#####
#####
# Script for plotting every solution
#####
#####
for i=1:length(results["files"])
    println(i)
    fig = plot_farm_figure(results["files"][i]);
    savefig("Farms_C32_figures/examples/solution_"results["files"][i][39:end-4]*".png");
    close(fig)
end

#####
#####
# Plot Histograms of the unit activation for each solution
#####
#####
examples,results = load("MiniC32_C32_examples_50_long.jld","examples","results");
function make_accuracy_histo(filename, examples, results, opto, proanti)
    i = find(results["files"] .== filename); 
    dex1 = vec(examples[i,opto,proanti,1,end,:]);
    dex2 = vec(examples[i,opto,proanti,4,end,:]);
    figure()
    subplot(2,1,1);
    plt[:hist](dex1,20,color=(0/255,128/255,0/255))
    ylim(0, 100)
    plot(vec([0, 1]), vec([0 0]), color=(179/255,179/255,179/255));
    subplot(2,1,2);
    plt[:hist](dex2,20,color=(100/255,200/255,100/255))
    ylim(0, 100)
    plot(vec([0, 1]), vec([0 0]), color=(179/255,179/255,179/255));
end



